#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 blast.grouped\n\n";

my $blast_out = $ARGV[0] or die $usage;

main: {


    my $counter = 0;

    my %query_to_top_hit; # only storing the hit with the greatest blast score.

    # outfmt6:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    
    open (my $fh, $blast_out) or die "Error, cannot open file $blast_out";
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $query_id = $x[0];
        my $db_id = $x[1];
        my $percent_id = $x[2];
        my $Evalue = $x[3];
        
        my $pct_target_len = $x[7];
        
        
        if ( (! exists $query_to_top_hit{$query_id}) || ($Evalue < $query_to_top_hit{$query_id}->{Evalue}) ) {
            
            $query_to_top_hit{$query_id} = { query_id => $query_id,
                                             db_id => $db_id,
                                             percent_id => $percent_id,
                                             Evalue => $Evalue,
                                             pct_target_len => $pct_target_len,
            };
            
        }
    }
    close $fh;
    
    
    ## get the best transcript hit per db ID
    
    my %db_id_to_trans_hits;
    
    foreach my $hit_struct (values %query_to_top_hit) {
        
        my $db_id = $hit_struct->{db_id};
        
        push (@{$db_id_to_trans_hits{$db_id}}, $hit_struct);
    }

        
    ## histogram summary

    my @bins = qw(10 20 30 40 50 60 70 80 90 100);
    my %bin_counts;

    
    foreach my $db_id (keys %db_id_to_trans_hits) {
        
        my @hit_structs = @{$db_id_to_trans_hits{$db_id}};
        
        @hit_structs = sort {$a->{Evalue} < $b->{Evalue}
                             ||
                                 $a->{pct_target_len} > $b->{pct_target_len} } @hit_structs;
        
        my $entry = shift @hit_structs; # take the lowest E-value w/ longest pct_target_len

        my $pct_cov = $entry->{pct_target_len};
        
        my $prev_bin = 0;
        foreach my $bin (@bins) {
            if ($pct_cov > $prev_bin && $pct_cov <= $bin) {
                $bin_counts{$bin}++;
            }
            $prev_bin = $bin;
        }
                
    }
   

    ## Report counts per bin
    print "#hit_pct_cov_bin\tcount_in_bin\t>bin_below\n";

    my $cumul = 0;
    foreach my $bin (reverse(@bins)) {
        my $count = $bin_counts{$bin} || 0;
        $cumul += $count;
        print join("\t", $bin, $count, $cumul) . "\n";
    }

    
    
    exit(0);
    

}
