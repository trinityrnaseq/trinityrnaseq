#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $help_flag;


my $MIN_PCT_MAX_EXPR = 5;

my $CPU = 2;

my $REF_MIN_PCT_LEN = 90;

my $usage = <<__EOUSAGE__;


Algorithm is as follows:

    Run CDHIT to cluster according to sequence identity and length

    A reference sequence is selected as having the highest expression level and within ${REF_MIN_PCT_LEN} % length  of the longest sequence in that cluster.

    Any sequences having less than $MIN_PCT_MAX_EXPR of the selected reference sequence are filtered out.

    Use the --filter_antisense_only flag to only have filtering applied to sequences with opposite transcribed orientation from the selected reference sequence in each cluster.

    Output is a new fasta file containing only those entries that pass the filtering criteria.

    
################################################################
#
#  --transcripts_fasta <string>        target transcript fasta file
#
#  --expr_matrix <string>              transcript expression matrix (TPM matrix)
#
#
# Optional:
#
#  --min_pct_max_expr <int>            default: $MIN_PCT_MAX_EXPR
#
#  --filter_antisense_only             default: off
#
#  --CPU <int>                         number of threads
#
#  --ref_min_pct_len <int>             minimum percent of max length to define candidate reference sequences
#                                      default: $REF_MIN_PCT_LEN
#
###########################################################################


__EOUSAGE__

    ;


# cd-hit-est  -i axo.indropB.iworm.fa -o axo.indropB.iworm.fa.cdhit98.fa -c 0.98 -aS 0.95 -T 20

my $transcripts_fasta_file;
my $expr_matrix_file;
my $filter_antisense_only;

my $DEBUG = 0;


&GetOptions ( 'h' => \$help_flag,
              'transcripts_fasta=s' => \$transcripts_fasta_file,
              'expr_matrix=s' => \$expr_matrix_file,
              'filter_antisense_only' => \$filter_antisense_only,

              'CPU=i' => \$CPU,
              'min_pct_max_expr=i' => \$MIN_PCT_MAX_EXPR,
              'REF_MIN_PCT_LEN=i' => \$REF_MIN_PCT_LEN,
    

              'd' => \$DEBUG,
    );


if ($help_flag) {
    die $usage;
}


unless($transcripts_fasta_file && $expr_matrix_file) {
    die $usage;
}


main: {

    my $cdhit_prefix = "$transcripts_fasta_file.cdhit";
    
    my $cmd = "cd-hit-est -i $transcripts_fasta_file -o $cdhit_prefix -d 0 -c 0.98 -aS 0.95 -T $CPU";
    unless (-e "$cdhit_prefix.ok") {
        &process_cmd($cmd);

        &process_cmd("touch $cdhit_prefix.ok");
    }

    my %top_expr_val = &get_top_expr_val($expr_matrix_file);
    
    my %accs_retain = &filter_noisy_transcripts("$cdhit_prefix.clstr", \%top_expr_val, $REF_MIN_PCT_LEN, 
                                                $MIN_PCT_MAX_EXPR, $filter_antisense_only);

    # output refined fasta file:
    

    my $fasta_reader = new Fasta_reader($transcripts_fasta_file);

    my $count_kept = 0;
    my $count_excluded = 0;
    
    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();

        if ($accs_retain{$acc}) {
            print $seq_obj->get_FASTA_format();
            delete $accs_retain{$acc};
            $count_kept++;
        }
        else {
            $count_excluded++;
        }
    }
    
    if (%accs_retain) {
        die "Error, didnt retrieve sequences for accession: " . Dumper(%accs_retain);
    }
    else {
        print STDERR "Done. Excluded $count_excluded  = " . sprintf("%.2f", $count_excluded / ($count_kept + $count_excluded) * 100) . "% of sequences\n";
    }
    
    exit(0);
}

####
sub get_top_expr_val {
    my ($expr_matrix_file) = @_;

    my %top_expr_val;
    
    open(my $fh, $expr_matrix_file) or die "Error, cannot open file: $expr_matrix_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $trans_acc = shift @x;
        @x = sort {$a<=>$b} @x;
        my $max_expr = pop @x;
        $top_expr_val{$trans_acc} = $max_expr;

    }

    close $fh;

    return(%top_expr_val);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret";
    }

    return;
}

####
sub filter_noisy_transcripts {
    my ($cdhit_clstr_file, $top_expr_vals_href, $ref_min_pct_len, $min_pct_max_expr, $filter_antisense_only) = @_;


    
    
    my %accs_want;

    my $select_clusters_want_sref = sub {
        
        my @structs = @_;
        

    
        @structs = reverse sort {$a->{len}<=>$b->{len}} @structs;

        my $longest_len = $structs[0]->{len};

        my @ref_structs;
        foreach my $struct (@structs) {
            if ($struct->{len} / $longest_len * 100 >= $ref_min_pct_len) {
                push (@ref_structs, $struct);
            }
        }

        # take the highest expr entry as official ref seq
        @ref_structs = reverse sort {$a->{expr}<=>$b->{expr}} @ref_structs;
        
        
        my $ref_seq_struct = shift @ref_structs;
        $accs_want{ $ref_seq_struct->{acc} } = 1;
        
        my $audit_text = "* ref selected as: $ref_seq_struct->{acc} "
            . " len: $ref_seq_struct->{len}"
            . " expr: $ref_seq_struct->{expr}"
            . " [$ref_seq_struct->{orient}]\n";
        

        print STDERR "* selected ref seq as: " . Dumper($ref_seq_struct) if $DEBUG;
        
    
        my $min_allowed_expr = $ref_seq_struct->{expr} * $min_pct_max_expr / 100;

        my $excluded_flag = 0;
        
        foreach my $struct (@structs) {
            
            if ($struct eq $ref_seq_struct) { next; }
        
            my $retain_flag = 0;
    
            if ($filter_antisense_only && $struct->{orient} eq $ref_seq_struct->{orient}) {
                # all good, not antisense
                $retain_flag = 1;
            }

            elsif ($struct->{expr} >= $min_allowed_expr) {
                # meets expression criteria
                $retain_flag = 1;
            }
            
            if ($retain_flag) {
                $accs_want{ $struct->{acc} } = 1;
                $audit_text .= "-keeping: $struct->{acc} len: $struct->{len} expr: $struct->{expr} [$struct->{orient}]\n";
            }
            else {
                print STDERR "-excluding " . Dumper($struct) if $DEBUG;
                $audit_text .= "-EXCLUDING: $struct->{acc} len: $struct->{len} expr: $struct->{expr} [$struct->{orient}]\n";
                
                $excluded_flag = 1;
            }
        }

        if ($excluded_flag) {
            print STDERR "$audit_text\n";
        }
        
    };

    
    my @cluster;
    open(my $fh, $cdhit_clstr_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^>/) {
            if (@cluster) {
                &$select_clusters_want_sref(@cluster);
            }
            @cluster = ();
        }
        else {
            my ($idx, $len, $trans_acc, $selected_star, $orient_pct) = split(/\s+/);
            $len =~ s/nt,//;
            $trans_acc =~ s/>//;
            $trans_acc =~ s/\.\.\.$//;

            my $expr = $top_expr_vals_href->{$trans_acc};
            unless (defined $expr) {
                die "Error, no expr value for $trans_acc";
            }
            
            my ($orient, $pct) = ('+', '.');
            if ($selected_star eq 'at') {
                ($orient, $pct) = split(/\//, $orient_pct);
            }
            
            my $struct = { acc => $trans_acc,
                           len => $len,
                           orient => $orient,
                           expr => $expr,
            };
            push (@cluster, $struct);
            
        }

    }
    close $fh;
    
    if (@cluster) {
        &$select_clusters_want_sref(@cluster);
    }

    return(%accs_want);
}
