#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Data::Dumper;


my $help_flag;


my $usage = <<__EOUSAGE__;

####################################################################################
#
# required:
#
#  --transcripts_fasta <string>   fasta file containing transcript sequences.      
#
#  --gene_trans_map <string>      gene-trans mapping file (format: gene_id(tab)transcript_id)
#
#  --expr_matrix <string>         expression matrix (fpkm or tpm, ideally TMM-normalized)
#
#  --samples_file <string>        samples file
#
#  --min_pct_iso <float>          minimum percent of total isoform expression
#
# optional:
#
#  --out_prefix <string>          output filename prefix ("pruned.\${min_pct_iso}"
#
####################################################################################

__EOUSAGE__

    ;



my $gene_trans_map_file;
my $expr_matrix;
my $samples_file;
my $min_pct_iso;
my $out_prefix;
my $transcripts_fasta_file;

&GetOptions ( 'h' => \$help_flag,

              'gene_trans_map=s' => \$gene_trans_map_file,
              'expr_matrix=s' => \$expr_matrix,
              'transcripts_fasta=s' => \$transcripts_fasta_file,

              'samples_file=s' => \$samples_file,
              'min_pct_iso=s' => \$min_pct_iso,
              'out_prefix=s' => \$out_prefix,
              
    );

unless ($gene_trans_map_file && $expr_matrix && $samples_file && $min_pct_iso && $transcripts_fasta_file) {
    die $usage;
}

unless ($out_prefix) {
    $out_prefix = "pruned.$min_pct_iso";
}

main: {

    print STDERR "-parsing samples file: $samples_file\n";
    my %sample_to_replicates = &parse_samples_file($samples_file);

    print STDERR "-parsing gene-trans mapping file: $gene_trans_map_file\n";
    my %gene_to_trans = &parse_gene_trans_relationships($gene_trans_map_file);
    
    print STDERR "-parsing expr matrix: $expr_matrix\n";
    my %trans_to_avg_expr = &parse_expr_matrix($expr_matrix, \%sample_to_replicates);
    
    my $out_log_filename = "$out_prefix.log";
    my $out_fasta_filename = "$out_prefix.fasta";

    open(my $out_log_fh, ">$out_log_filename") or die "Error, cannot write to $out_log_filename";
    open(my $out_fasta_fh, ">$out_fasta_filename") or die "Error, cannot write to $out_fasta_filename";

    print STDERR "-pruning isoforms\n";

    ## decide on which transcripts to keep
    my %trans_retain;
    foreach my $gene (keys %gene_to_trans) {
        my @trans_ids = keys %{$gene_to_trans{$gene}};
        if (scalar(@trans_ids) == 1) {
            # no alt splicing, just keep it.
            $trans_retain{ $trans_ids[0] } = 1;
            next;
        }

        my %trans_to_sum_expr_all_conditions;
        my @trans_to_sample_expr;
        my %sample_to_sum_expr;
        foreach my $sample_name (keys %sample_to_replicates) {
            my $sum_expr = 0;
            my @trans_to_expr;
            foreach my $trans_id (@trans_ids) {
                my $expr = $trans_to_avg_expr{$trans_id}->{$sample_name};
                unless (defined $expr) {
                    die "Error, no expression set for trans_id: [$trans_id] and sample type: [$sample_name] ";
                }
                $trans_to_sum_expr_all_conditions{$trans_id} += $expr;
                push (@trans_to_sample_expr, [$trans_id, $sample_name, $expr]);
                $sample_to_sum_expr{$sample_name} += $expr;
            }
        }

        # prioritize transcript by order across experiments
        @trans_to_sample_expr = reverse sort {$trans_to_sum_expr_all_conditions{$a->[0]} <=> $trans_to_sum_expr_all_conditions{$b->[0]}} @trans_to_sample_expr;
        
        my %seen_sample;
        foreach my $trans_info_aref (@trans_to_sample_expr) {
            # always retain top isoform in any sample type.
            my ($trans_id, $sample_name, $expr) = @$trans_info_aref;
            
            $expr = sprintf("%.3f", $expr);
            
            my $retain_flag = 0;
            my $selected_top_entry = 0;
            if (! $seen_sample{$sample_name}) {
                # top expr entry for that sample
                $seen_sample{$sample_name} = 1;
                $selected_top_entry = 1;
                $retain_flag = 1;
            }
            
            my $sum_sample_expr = $sample_to_sum_expr{$sample_name};
            my $pct_iso = 0;
            if ($sum_sample_expr > 0) {
                $pct_iso = sprintf("%.2f", $expr / $sum_sample_expr * 100);
                if ($pct_iso >= $min_pct_iso) {
                    $retain_flag = 1;
                }
            }
            
            if ($retain_flag) {
                $trans_retain{$trans_id} = 1;
            }
            
            ## output log entry
            my $out_text = join("\t", $gene, $trans_id, $sample_name, $expr, $pct_iso);
            if ($retain_flag) {
                $out_text .= "\t[RETAINED]";
            }
            if ($selected_top_entry) {
                $out_text .= "\t[TOP_EXPR]";
            }
            print $out_log_fh "$out_text\n";
        }
        
        print $out_log_fh "\n"; # spacer between genes.
    }
    

    close $out_log_fh;
    
    my $fasta_reader = new Fasta_reader($transcripts_fasta_file) or die $!;
    my %to_retrieve = %trans_retain;

    print STDERR "-outputting filtered transcripts as: $out_fasta_filename\n";
    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();
        if ($trans_retain{$acc}) {
            my $fasta_entry = $seq_obj->get_FASTA_format();
            chomp $fasta_entry;
            print $out_fasta_fh "$fasta_entry\n";
            delete($to_retrieve{$acc});
        }
    }
    
    if (%to_retrieve) {
        die "Failed to retrieve fasta records for transcripts: " . Dumper(\%to_retrieve);
    }
    
    close $out_fasta_fh;
    
    print STDERR "-done.  See outputs: $out_prefix.\*\n\n";
    
    exit(0);
}

####
sub parse_gene_trans_relationships {
    my ($gene_trans_map_file) = @_;

    my %gene_to_trans;

    open(my $fh, $gene_trans_map_file) or die "Error, cannot open file $gene_trans_map_file";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        chomp;
        
        my ($gene_id, $trans_id) = split(/\t/);
        $gene_to_trans{$gene_id}->{$trans_id} = 1;
    }
    close $fh;
    
    return(%gene_to_trans);

}


####
sub parse_samples_file {
    my ($filename) = @_;

    my %sample_to_reps;
    
    open (my $fh, $filename) or die "Error, cannot open file $filename";
    while (<$fh>) {
        chomp;
        my ($condition, $replicate_name, @rest) = split(/\t/);
        
        $sample_to_reps{$condition}->{$replicate_name} = 1;
    }

    return(%sample_to_reps);
}


####
sub parse_expr_matrix {
    my ($expr_matrix, $sample_to_replicates_href) = @_;

    open(my $fh, $expr_matrix) or die "Error, cannot open file $expr_matrix";
    
    my $header = <$fh>;
    chomp $header;
    my @header_fields = split(/\t/, $header);
    
    my %replicate_to_col_number;

    my %trans_expr;

    my $first_data_row = 1;
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);

        # process column info as needed.
        if ($first_data_row) {
            if (scalar(@header_fields) == scalar(@x) - 1) {
                # adjust header
                unshift(@header_fields, '');
            }
            # assign header column numbers
            for (my $i = 1; $i <= $#header_fields; $i++) {
                my $rep_name = $header_fields[$i];
                $replicate_to_col_number{$rep_name} = $i;
            }
            $first_data_row = 0;
        }
        
        my $trans_id = $x[0];
        
        ## process data, compute averages
        foreach my $sample_name (keys %$sample_to_replicates_href) {
            my @expr_vals;
            my @rep_names = keys %{$sample_to_replicates_href->{$sample_name}};
            foreach my $rep_name (@rep_names) {
                my $col_id = $replicate_to_col_number{$rep_name} || die "Error, no column identified for replicate name: [$rep_name] ";
                my $expr_val = $x[$col_id];
                push (@expr_vals, $expr_val);
            }
            my $avg_expr = &sum(@expr_vals) / scalar(@expr_vals);

            $trans_expr{$trans_id}->{$sample_name} = $avg_expr;
        }
    }

    return(%trans_expr);
}

####
sub sum {
    my @vals = @_;

    my $sum_val = 0;
    foreach my $v (@vals) {
        $sum_val += $v;
    }
    
    return($sum_val);
}
