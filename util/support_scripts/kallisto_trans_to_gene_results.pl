#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "\n\nusage: $0 abundance.tsv  gene_to_trans_map_file.txt\n\n\n";

my $abundance_tsv = $ARGV[0] or die $usage;
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
    
    
    open (my $fh, $abundance_tsv) or die "Error, cannot open file $abundance_tsv";
    my $header = <$fh>;
    chomp $header;
    my %field_index;
    my @fields = split(/\t/, $header);
    {
        
        for (my $i = 0; $i <= $#fields; $i++) {
            my $field = $fields[$i];
            $field_index{$field} = $i;
        }
    }
    

    my %gene_data;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        my $trans_id = $x[ $field_index{target_id} ];
        my $tpm = $x[ $field_index{tpm} ];
        my $length = $x[ $field_index{length} ];
        my $eff_length = $x[ $field_index{eff_length} ];
        my $est_counts = $x[ $field_index{est_counts} ];
        
        my $gene = $trans_to_gene_info{$trans_id} or die "Error, cannot find gene identifier for transcript [$trans_id] ";
        
        push (@{$gene_data{$gene}}, { trans_id => $trans_id, 
                                      tpm => $tpm,
                                      length => $length,
                                      eff_length => $eff_length,
                                      est_counts => $est_counts,
                                  });
        

    }
    close $fh;


    ## Output gene summaries:
    
    print $header . "\n";

    foreach my $gene (keys %gene_data) {
        my @trans_structs = @{$gene_data{$gene}};

        my @trans_ids;
        my $sum_counts = 0;
        my $sum_tpm = 0;
        
        my $counts_per_len_sum = 0;
        my $counts_per_eff_len_sum = 0;
        
        my $sum_lengths = 0;
        my $sum_eff_lengths = 0;
        
        my $num_trans = scalar(@trans_structs);


        foreach my $struct (@trans_structs) {
            
            #print Dumper($struct);

            my $trans_id = $struct->{trans_id};
            my $tpm = $struct->{tpm};
            my $length = $struct->{length};
            
            my $eff_length = $struct->{eff_length};
            my $est_counts = $struct->{est_counts};
            
            unless ($eff_length > 0) {
                $eff_length = 1; # cannot have zero length feature!
            }
            
            unless ($length > 0 && $eff_length > 0) {
                die "Error, length: $length, eff_length: $eff_length" . Dumper($struct);
            }

            $sum_lengths += $length;
            $sum_eff_lengths += $eff_length;

            $counts_per_len_sum += $est_counts/$length;
            
            $counts_per_eff_len_sum += $est_counts/$eff_length;
            
            $sum_counts += $est_counts;
            $sum_tpm += $tpm;
        }
        
        my $gene_length = $sum_lengths / $num_trans;
        my $gene_eff_length = $sum_eff_lengths / $num_trans;
        if ($sum_counts) {
            # set lengths as weighted by expression of isoforms.
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
                          tpm => sprintf("%.2f", $sum_tpm),
                          length => sprintf("%.2f", $gene_length),
                          eff_length => sprintf("%.2f", $gene_eff_length),
                          est_counts => sprintf("%.2f", $sum_counts),
                          
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
