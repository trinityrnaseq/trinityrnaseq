#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "\n\nusage: $0 results.xprs  gene_to_trans_map_file.txt\n\n\n";

my $results_xprs = $ARGV[0] or die $usage;
my $gene_to_trans_map_file = $ARGV[1] or die $usage;


main: {
    
    my %trans_to_gene_info;
    {
        open (my $fh, $gene_to_trans_map_file) or die "Error, cannot open file $gene_to_trans_map_file";
        while (<$fh>) {
            unless (/\w/) { next; }
            chomp;
            my ($gene, $trans, @rest) = split(/\s+/);
            unless ($gene && $trans) {
                die "Error, cannot extract gene & trans relationship from line $_ of file $gene_to_trans_map_file";
            }
            $trans_to_gene_info{$trans} = $gene;
        }
        close $fh;
    }
    
    
    open (my $fh, $results_xprs) or die "Error, cannot open file $results_xprs";
    my $header = <$fh>;
    chomp $header;
    my %field_index;
    my @fields = split(/\t/, $header);
    {
        
        for (my $i = 1; $i <= $#fields; $i++) {
            my $field = $fields[$i];
            $field_index{$field} = $i;
        }
    }
    

    my %gene_data;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        my $trans_id = $x[ $field_index{target_id} ];
        my $fpkm = $x[ $field_index{fpkm} ];
        my $length = $x[ $field_index{length} ];
        my $eff_length = $x[ $field_index{eff_length} ];
        my $eff_counts = $x[ $field_index{eff_counts} ];
        my $tpm = $x[ $field_index{tpm} ];
    
        my $gene = $trans_to_gene_info{$trans_id} or die "Error, cannot find gene identifier for transcript [$trans_id] ";
        
        push (@{$gene_data{$gene}}, { trans_id => $trans_id, 
                                      fpkm => $fpkm,
                                      length => $length,
                                      eff_length => $eff_length,
                                      eff_counts => $eff_counts,
                                      tpm => $tpm,
              });
        

    }
    close $fh;


    ## Output gene summaries:
    
    print $header . "\n";

    foreach my $gene (keys %gene_data) {
        my @trans_structs = @{$gene_data{$gene}};

        my @trans_ids;
        my $sum_counts = 0;
        my $sum_fpkm = 0;
        my $sum_tpm = 0;

        my $counts_per_len_sum = 0;
        my $counts_per_eff_len_sum = 0;
        
        my $sum_lengths = 0;
        my $sum_eff_lengths = 0;
        
        my $num_trans = scalar(@trans_structs);
        
        foreach my $struct (@trans_structs) {
            
            #print Dumper($struct);

            my $trans_id = $struct->{trans_id};
            my $fpkm = $struct->{fpkm};
            my $tpm = $struct->{tpm};
            my $length = $struct->{length};
            

            my $eff_length = $struct->{eff_length};
            my $eff_counts = $struct->{eff_counts};
            
            unless ($eff_length > 0) {
                $eff_length = 1; # cannot have zero length feature!
            }
            
            unless ($length > 0 && $eff_length > 0) {
                die "Error, length: $length, eff_length: $eff_length" . Dumper($struct);
            }

            $sum_lengths += $length;
            $sum_eff_lengths += $eff_length;

            $counts_per_len_sum += $eff_counts/$length;
            
            $counts_per_eff_len_sum += $eff_counts/$eff_length;
            
            $sum_counts += $eff_counts;
            $sum_fpkm += $fpkm;

            $sum_tpm += $tpm;
        }
        
        my $gene_length = $sum_lengths / $num_trans;
        my $gene_eff_length = $sum_eff_lengths / $num_trans;
        if ($sum_counts) {
            # compute it based on the formula:  FPKM_gene = FPKM_isoA + FPKM_IsoB + ... 
            eval {
                $gene_length = $sum_counts / $counts_per_len_sum;
                $gene_eff_length = $sum_counts / $counts_per_eff_len_sum;
            };
            if ($@) {
                print STDERR "$@\n" . Dumper(\@trans_structs);
                die;
            }
        }
        
        my %gene_info = ( target_id => $gene,
                          fpkm => sprintf("%.2f", $sum_fpkm),
                          length => sprintf("%.2f", $gene_length),
                          eff_length => sprintf("%.2f", $gene_eff_length),
                          eff_counts => sprintf("%.2f", $sum_counts),
                          tpm => $sum_tpm,
                          );
        
        my @vals;
        foreach my $field (@fields) {
            my $result = $gene_info{$field};
            unless (defined $result) {
                $result = "NA";
            }
            push (@vals, $result);
        }
        print join("\t", @vals) . "\n";
    }
    

    exit(0);
}
