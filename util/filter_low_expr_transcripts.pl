#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use Fasta_reader;

my $help_flag;


my $usage = <<__EOUSAGE__;

##########################################################################################
#
#  --matrix|m <string>            expression matrix (TPM or FPKM, *not* raw counts)
#
#  --transcripts|t <string>       transcripts fasta file (eg. Trinity.fasta)
#
#
#  # expression level filter:
#
#     --min_expr_any <float>      minimum expression level required across any sample (default: 0)
#
#  # Isoform-level filtering
#
#     --min_pct_dom_iso <int>         minimum percent of dominant isoform expression (default: 0)
#          or
#     --highest_iso_only          only retain the most highly expressed isoform per gene (default: off)
#                                 (mutually exclusive with --min_pct_dom_iso param)
#
#     # requires gene-to-transcript mappings
#
#     --trinity_mode              targets are Trinity-assembled transcripts
#         or
#     --gene_to_trans_map <string>   file containing gene-to-transcript mappings
#                                    (format is:   gene(tab)transcript )
#
#########################################################################################


__EOUSAGE__

    ;


my $matrix_file;
my $transcripts_file;
my $min_expr_any = 0;
my $min_pct_dom_iso = 0;
my $highest_iso_only_flag = 0;
my $trinity_mode_flag = 0;
my $gene_to_trans_map_file;


&GetOptions ( 'help|h' => \$help_flag,

              'matrix|m=s' => \$matrix_file,
              'transcripts|t=s' => \$transcripts_file,
              
              'min_expr_any=f' => \$min_expr_any,
              'min_pct_dom_iso=i' => \$min_pct_dom_iso,
              'highest_iso_only' => \$highest_iso_only_flag,

              'trinity_mode' => \$trinity_mode_flag,
              'gene_to_trans_map=s' => \$gene_to_trans_map_file,
              

    );


if ($help_flag) {
    die $usage;
}


unless ($matrix_file && $transcripts_file &&
        ($min_expr_any || $min_pct_dom_iso || $highest_iso_only_flag) ) {

    die $usage;
}

if ( ($min_pct_dom_iso || $highest_iso_only_flag) && ! ($trinity_mode_flag || $gene_to_trans_map_file) ) {
    die "Error, if --min_pct_dom_iso or --highest_iso_only, must also specify either --trinity_mode or --gene_to_trans_map";
}

if ($min_pct_dom_iso && $highest_iso_only_flag) {
    die "Error, --min_pct_dom_iso and --highest_iso_only are mutually exclusive parameters. ";
}


main: {

    my %expr_vals = &parse_expr_matrix($matrix_file);
    
    if ($min_pct_dom_iso || $highest_iso_only_flag) {
        
        my %gene_to_iso_map = ($trinity_mode_flag) 
            ? &parse_Trinity_gene_mapping($transcripts_file)
            : &parse_gene_trans_map_file($gene_to_trans_map_file);
        
        &add_pct_iso_stats(\%expr_vals, \%gene_to_iso_map);
    }

    my $total_records = 0;
    my $retained_records = 0;
    
    my $fasta_reader = new Fasta_reader($transcripts_file);
    while (my $seq_obj = $fasta_reader->next()) {

        $total_records++;
        
        my $acc = $seq_obj->get_accession();

        my $keep_flag = 1;

        my $info_struct = $expr_vals{$acc} or die "Error, no expression record stored for acc: [$acc].  Be sure to provide the transcript expression matrix and all transcripts in the $transcripts_file must have records in the transcript expression matrix file.";
        
        if ($min_expr_any && $info_struct->{max_expr} < $min_expr_any) {
            $keep_flag = 0;
            print STDERR "-excluding $acc, max_expr: $info_struct->{max_expr} < $min_expr_any\n";
        }
        if ($min_pct_dom_iso && (! $info_struct->{top_iso}) && $info_struct->{pct_dom_iso_expr} < $min_pct_dom_iso) {
            # notice we'll still keep the dominant isoform for the gene even if it's pct iso < $min_pct_dom_iso.
            ## dont want to be silly and throw out the gene altogther...  :)
            print STDERR "-excluding $acc, pct_dom_iso_expr $info_struct->{pct_dom_iso_expr} < $min_pct_dom_iso\n";
            
            $keep_flag = 0;
        }
        
        if ($highest_iso_only_flag && ! $info_struct->{top_iso}) {
            print STDERR "-excluding $acc, not top_iso\n";
            $keep_flag = 0;
        }

        if ($keep_flag) {
            $retained_records++;
            my $fasta_record = $seq_obj->get_FASTA_format();
            chomp $fasta_record;
            my ($header_line, @seq_lines) = split(/\n/, $fasta_record);
            # tack on the pct expr info onto the header
            my $top_iso_flag = $info_struct->{top_iso};
            my $pct_iso_expr = (defined $info_struct->{pct_iso_expr}) ? $info_struct->{pct_iso_expr} : "NA";
            
            my $pct_dom_iso_expr = (defined $info_struct->{pct_dom_iso_expr}) ? $info_struct->{pct_dom_iso_expr} : "NA";
            
            $header_line .= " top_iso:$top_iso_flag pct_iso_expr=$pct_iso_expr pct_dom_iso_expr=$pct_dom_iso_expr max_expr_any=$info_struct->{max_expr}";
            
            print join("\n", $header_line, @seq_lines) . "\n";
            
        }
    }
    
    my $pct_records_retained = sprintf("%.2f", $retained_records / $total_records * 100);
    print STDERR "\n\n\tRetained $retained_records / $total_records = $pct_records_retained\% of total transcripts.\n\n\n";
    
    
    exit(0);
    
    
}

####
sub add_pct_iso_stats {
    my ($expr_vals_href, $gene_to_iso_map_href) = @_;
    
    foreach my $gene (keys %$gene_to_iso_map_href) {
        
        my @isoforms = keys %{$gene_to_iso_map_href->{$gene}};
        
        if (scalar @isoforms == 1) {
            # only one isoform, so must be 100% of that gene.
            $expr_vals_href->{ $isoforms[0] }->{pct_iso_expr} = 100;
            $expr_vals_href->{ $isoforms[0] }->{pct_dom_iso_expr} = 100;
            $expr_vals_href->{ $isoforms[0] }->{top_iso} = 1;
            
        }
        else {
            # determine fraction of total gene expr
            # first, get sum of gene expr across isoforms
            my $gene_sum_expr = 0;
            my $dominant_iso_expr = 0;
            foreach my $iso (@isoforms) {
                
                my $expr = $expr_vals_href->{$iso}->{sum_expr};
                if (!defined($expr)) {
                    use Data::Dumper;
                    print STDERR "ISO: $iso\t" . Dumper($expr_vals_href->{$iso});
                }
                if ($expr > $dominant_iso_expr) {
                    $dominant_iso_expr = $expr;
                }

                
                $gene_sum_expr += $expr;
            }
            # now compute pct iso
            foreach my $iso (@isoforms) {
                my $expr = $expr_vals_href->{$iso}->{sum_expr};
                my $pct_iso = 0;
                if ($gene_sum_expr > 0) {
                    $pct_iso = sprintf("%.2f", $expr / $gene_sum_expr * 100);
                }
                $expr_vals_href->{$iso}->{pct_iso_expr} = $pct_iso;

                my $pct_dom_iso_expr = 0;
                if ($dominant_iso_expr > 0) {
                    $pct_dom_iso_expr = sprintf("%.2f", $expr / $dominant_iso_expr * 100);
                }
                $expr_vals_href->{$iso}->{pct_dom_iso_expr} = $pct_dom_iso_expr;
            }
            # set top iso
            @isoforms = sort { $expr_vals_href->{$a}->{pct_iso_expr} <=> $expr_vals_href->{$b}->{pct_iso_expr} } @isoforms;
            
            my $top_isoform = pop @isoforms;  # note, if there's no gene expression for some reason, choice isn't informative.
            $expr_vals_href->{$top_isoform}->{top_iso} = 1;
        }
    }
}


####
sub parse_expr_matrix {
    my ($matrix_file) = @_;
    
    my %expr_vals;
    
    open (my $fh, $matrix_file) or die "Error, cannot open file $matrix_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @expr = split(/\t/);
        my $acc = shift @expr;

        my $max_val = 0;
        my $sum = 0;
        
        foreach my $expr_val (@expr) {
            $sum += $expr_val;
            if ($expr_val > $max_val) {
                $max_val = $expr_val;
            }
        }

        $expr_vals{$acc}->{max_expr} = $max_val;
        $expr_vals{$acc}->{sum_expr} = $sum;
        $expr_vals{$acc}->{pct_iso_expr} = undef; # set later
        $expr_vals{$acc}->{pct_dom_iso_expr} = undef;
        $expr_vals{$acc}->{top_iso} = 0; # set later to the isoform with highest expression for that gene.
        
    }
    close $fh;
    
    return(%expr_vals);
}

####
sub parse_Trinity_gene_mapping {
    my ($transcripts_file) = @_;

    my %gene_to_iso_map;
    
    open (my $fh, $transcripts_file) or die "Error, cannot open file $transcripts_file";
    while (<$fh>) {
        if (/^>(\S+)/) {
            my $acc = $1;
            $acc =~ /^(\S+)(_i\d+)$/ or die "Error, cannot parse Trinity accession: $acc";
            my $gene_id = $1;
            
            $gene_to_iso_map{$gene_id}->{$acc} = 1;
        }
    }
    close $fh;

    return(%gene_to_iso_map);
}

####
sub parse_gene_trans_map_file {
    my ($gene_to_trans_map_file) = @_;

    my %gene_to_iso_map;
    
    open (my $fh, $gene_to_trans_map_file) or die "Error, cannot open file $gene_to_trans_map_file";
    while (<$fh>) {
        chomp;
        my ($gene, $trans) = split(/\t/);
        
        $gene_to_iso_map{$gene}->{$trans} = 1;
    }
    close $fh;

    return (%gene_to_iso_map);
}

