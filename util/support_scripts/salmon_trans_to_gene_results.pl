#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "\n\nusage: $0 quant.sf gene_to_trans_map_file.txt\n\n\n";

my $quant_sf = $ARGV[0] or die $usage;
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
    
    
    open (my $fh, $quant_sf) or die "Error, cannot open file $quant_sf";
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
        
        # quant.sf format:
        #
        #Name    Length  EffectiveLength TPM     NumReads
        #TRINITY_DN10_c0_g1_i1   334     67.2849 3125.31 7
        #TRINITY_DN11_c0_g1_i1   319     55.1277 0       0
        #TRINITY_DN12_c0_g1_i1   244     244     1231.18 10
        #TRINITY_DN17_c0_g1_i1   229     229     393.549 3
        #TRINITY_DN18_c0_g1_i1   633     360.371 593.619 7.12107
        
        my @x = split(/\t/);
        
        my $trans_id = $x[ $field_index{Name} ];
        my $tpm = $x[ $field_index{TPM} ];
        my $length = $x[ $field_index{Length} ];
        my $eff_length = $x[ $field_index{EffectiveLength} ];
        my $est_counts = $x[ $field_index{NumReads} ];
        
        my $gene = $trans_to_gene_info{$trans_id} or die "Error, cannot find gene identifier for transcript [$trans_id] ";
        
        push (@{$gene_data{$gene}}, { Name => $trans_id, 
                                      TPM => $tpm,
                                      Length => $length,
                                      EffectiveLength => $eff_length,
                                      NumReads => $est_counts,
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

            my $trans_id = $struct->{Name};
            my $tpm = $struct->{TPM};
            my $length = $struct->{Length};
            
            my $eff_length = $struct->{EffectiveLength};
            my $est_counts = $struct->{NumReads};
            
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
        
        my %gene_info = ( Name => $gene,
                          TPM => sprintf("%.2f", $sum_tpm),
                          Length => sprintf("%.2f", $gene_length),
                          EffectiveLength => sprintf("%.2f", $gene_eff_length),
                          NumReads => sprintf("%.2f", $sum_counts),
                          
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
