#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);



my $help_flag;
my $FDR_thresh = 0.05;


my $usage = <<__EOUSAGE__;

###############################################################
#
#
#  Required:
#
#  --trans_DE_results_dir <string>   /path/to/DE_results/dir/  
#                                      (contains the *.DE_results files) from pairwise comparisons.
#
#  --trinity_fasta <string>          Trinity fasta file
#
#  Optional:
#
#  --max_FDR <float>                 maximum FDR (default: 0.05)
#
##############################################################


__EOUSAGE__

    ;



my $trans_results_dir = "";
my $trinity_fasta_file = "";

&GetOptions ( 'h' => \$help_flag, 
              'trans_DE_results_dir=s' => \$trans_results_dir,
              'trinity_fasta=s' => \$trinity_fasta_file,
              
              'max_FDR=f' => \$FDR_thresh,
              );


if ($help_flag) {
    die $usage;
}

unless ($trans_results_dir && $trinity_fasta_file) {
    die $usage;
}



# want DE isoforms
# do we see isoform switching?
# multiple isoforms found as DE between each pair of conditions, and opposite directions

my @iso_swap_pairs;


main: {

    my %isoforms_per_gene = &parse_isoforms_per_gene($trinity_fasta_file);
    
    
    my %cond_pair_to_trans_results = &parse_DE_result_files($trans_results_dir);
    
    print join("\t", "iso_per_gene", "cond_pair", "score", 
               "feature", "logFC", "logCPM", "FDR",
               "feature", "logFC", "logCPM", "FDR") . "\n";
    
    
    foreach my $cond_pair (keys %cond_pair_to_trans_results) {
        my @DE_multi_iso = &get_DE_genes($cond_pair_to_trans_results{$cond_pair}, $FDR_thresh);
        
            
        foreach my $iso_set (@DE_multi_iso) {
            my @iso_structs = @$iso_set;
            
            my @iso_meet_fdr_thresh;
            foreach my $iso_struct (@iso_structs) {
                if ($iso_struct->{FDR} <= $FDR_thresh) {
                    push (@iso_meet_fdr_thresh, $iso_struct);
                }
            }
            
            if (scalar @iso_meet_fdr_thresh > 1) {
                &examine_for_iso_swap($cond_pair, $FDR_thresh, \@iso_meet_fdr_thresh);
            }
        }
    }
    

    
    @iso_swap_pairs = reverse sort {$a->{score}<=>$b->{score}} @iso_swap_pairs;

    foreach my $iso_swap (@iso_swap_pairs) {

        my $trans_name = $iso_swap->{A}->{feature};
        my $iso_count_for_gene;
        if ($trans_name =~ /^(.*c\d+_g\d+)_i\d+/) {
            my $gene_id = $1;
            $iso_count_for_gene = $isoforms_per_gene{$gene_id};
        }
        else {
            die "Error, cannot extract gene info from $trans_name";
        }
        
        print join("\t", $iso_count_for_gene, $iso_swap->{cond_pair}, sprintf("%.2f", $iso_swap->{score}),
                   $iso_swap->{A}->{feature}, $iso_swap->{A}->{logFC}, $iso_swap->{A}->{logCPM}, $iso_swap->{A}->{FDR},
                   $iso_swap->{B}->{feature}, $iso_swap->{B}->{logFC}, $iso_swap->{B}->{logCPM}, $iso_swap->{B}->{FDR},
                   ) . "\n";
    }
    

    exit(0);

}

use Data::Dumper;

####
sub examine_for_iso_swap {
    my ($cond_pair, $FDR_thresh, $iso_set_aref) = @_;

    for (my $i = 0; $i < $#$iso_set_aref; $i++) {

        for (my $j = $i + 1; $j <= $#$iso_set_aref; $j++) {

            my ($logFC_A, $log_FC_B) = sort {$a<=>$b} ($iso_set_aref->[$i]->{logFC}, $iso_set_aref->[$j]->{logFC});
            
            if ($logFC_A < 0 && $log_FC_B > 0) { 
                ## opposite polarity of change
                
                my $score = (-1*log($iso_set_aref->[$i]->{FDR})  + -1*log($iso_set_aref->[$j]->{FDR})) /2;

                push (@iso_swap_pairs,  { A => $iso_set_aref->[$i],
                                          B => $iso_set_aref->[$j],
                                          
                                          cond_pair => $cond_pair,
                                          
                                          score => $score,
                                      });

            }
        }
    }
}
                




####
sub get_DE_genes {
    my ($structs_aref, $FDR_thresh) = @_;

    my %genes_to_iso_structs;

    foreach my $struct (@$structs_aref) {

        if ($struct->{FDR} <= $FDR_thresh && abs($struct->{logFC}) > 2) {
            my $id = $struct->{feature};

            $id =~ /^(c\d+_g\d+)/ or die "Error, cannot extract trinity gene id from $id";
            
            my $gene_id = $1;
            $struct->{gene_id} = $gene_id;
            
            push (@{$genes_to_iso_structs{$gene_id}}, $struct);
                        
        }
    }
    
    
    my @multi_iso_DE_genes;

    foreach my $iso_set (values %genes_to_iso_structs) {
        
        if (scalar @$iso_set > 1) {
            push (@multi_iso_DE_genes, $iso_set);
        }
    }

    return(@multi_iso_DE_genes);
    

}

####
sub parse_DE_result_files {
    my ($dir) = @_;

    my %results;

    my @files = <$dir/*.DE_results>;

    foreach my $file (@files) {
    
        if ($file =~ /gill|embr/i) { next; }
        
        print STDERR "-processing file: $file\n";
        
        $file =~ /(\w+_vs_\w+)/ or die "Error, cannot extract pair of conditions from file $file";
        
        my $cond_pair = $1;

        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>;
        while (<$fh>) {
            chomp;
            my ($feature, $logFC, $logCPM, $PValue, $FDR) = split(/\t/);
            if ($FDR <= 0.05) {
                
                push (@{$results{$cond_pair}}, { feature => $feature,
                                                 logFC => $logFC,
                                                 logCPM => $logCPM,
                                                 PValue => $PValue,
                                                 FDR => $FDR,
                                             });
            }
        }
        
        close $fh;


    }


    return(%results);
}

####
sub parse_isoforms_per_gene {
    my ($trinity_fasta_file) = @_;

    my %iso_counter;

    open (my $fh, $trinity_fasta_file) or die $!;
    while (<$fh>) {
        if (/>(.*c\d+_g\d+)_i\d+/) {
            my $gene = $1;
            
            $iso_counter{$gene}++;
        }
    }
    close $fh;

    return(%iso_counter);
}

