#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib");
use Fasta_reader;
use Data::Dumper;

=ExampleCommands

# make blastable
makeblastdb \
  -in refTranscripts.fasta\
      -out refTranscripts -dbtype nucl

# run blast+
blastn -query Trinity.fasta -db refTranscripts -out blastn.fmt6.txt \
    -evalue 1e-20 -dust no -task megablast -num_threads 2 -max_target_seqs 1 -outfmt 6

# analyze results
analyze_blastPlus_topHit_coverage.pl blastn.fmt6.txt refTranscripts.fasta Trinity.fasta

=cut

    ;


my $usage = "usage: $0 blast+.outfmt6.txt query.fasta search_db.fasta [output_prefix=NameOfBlastFileHere] [verbose=0]\n\n";

my $blast_out = $ARGV[0] or die $usage;
my $fasta_file_A = $ARGV[1] or die $usage;
my $fasta_file_B = $ARGV[2] or die $usage; # the fasta files don't have to be in any special order.
my $output_prefix = $ARGV[3] || "$blast_out";
my $verbose = $ARGV[4] || 0;

main: {


    my $counter = 0;

    my %query_to_top_hit; # only storing the hit with the greatest blast score.

    # outfmt6:
    # qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    
    print STDERR "-parsing blast output: $blast_out\n" if $verbose;
    open (my $fh, $blast_out) or die "Error, cannot open file $blast_out";
    while (<$fh>) {
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $query_id = $x[0];
        my $db_id = $x[1];
        my $percent_id = $x[2];
        my $query_start = $x[6];
        my $query_end = $x[7];
        my $db_start = $x[8];
        my $db_end = $x[9];

        my $evalue = $x[10];
        my $bitscore = $x[11];
        
        if ( (! exists $query_to_top_hit{$query_id}) || ($bitscore > $query_to_top_hit{$query_id}->{bitscore}) ) {
            
            $query_to_top_hit{$query_id} = { query_id => $query_id,
                                             db_id => $db_id,
                                             percent_id => $percent_id,
                                             query_start => $query_start,
                                             query_end => $query_end,
                                             db_start => $db_start,
                                             db_end => $db_end,
                                             evalue => $evalue,
                                             bitscore => $bitscore,
                                         
                                             query_match_len => abs($query_end - $query_start) + 1,
                                             db_match_len => abs($db_end - $db_start) + 1,
                                             
                                             line => $line,
                                         };

        }

        $counter++;
        if ($counter % 100 == 0) {
            print STDERR "\r[$counter]    " if $verbose;
        }


    }
    close $fh;
    $counter = 0;
    print STDERR "\n" if $verbose;

    ## identify those entries we need sequence length info for.
    my %seq_lengths;
    my %seq_headers;
    {
        foreach my $entry (values %query_to_top_hit) {
            my $query_id = $entry->{query_id};
            my $db_id = $entry->{db_id};
            
            $seq_lengths{$query_id} = undef;
            $seq_lengths{$db_id} = undef;
        }
        
        ## get sequence length info
        foreach my $fasta_file ($fasta_file_A, $fasta_file_B) {
            
            print STDERR "-parsing seq length info from file: $fasta_file\n" if $verbose;
            
            my $fasta_reader = new Fasta_reader($fasta_file);
            
            while (my $seq_obj = $fasta_reader->next()) {
                
                $counter++;
                if ($counter % 100 == 0) {
                    print STDERR "\r[$counter]   " if $verbose;
                }


                my $acc = $seq_obj->get_accession();
                if (exists $seq_lengths{$acc}) {
                    
                    my $sequence = $seq_obj->get_sequence();
                    $seq_lengths{$acc} = length($sequence);
                
                    my $header = $seq_obj->get_header();
                    # remove the accession
                    my @header_pieces = split(/\s+/, $header);
                    shift @header_pieces;
                    $header = join(" ", @header_pieces);
                    $seq_headers{$acc} = $header;
                }
            }
            
            $counter = 0;
            print STDERR "\n" if $verbose;
                
        }
        
    }
    

    ## analyze the results.
    ## make this hit-centric, only retain the longest-coverage query hit for each database sequence.
    
    print STDERR "-analyzing hits.\n" if $verbose;
    my %db_id_to_greatest_pct_cov; # ties broken by bitscore
    {
        

        open (my $ofh, ">$output_prefix.w_pct_hit_length") or die $!;
        
        print $ofh join("\t", "#qseqid", "sseqid", "pident", "length", "mismatch",
                        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                        "db_hit_len", "pct_hit_len_aligned", "hit_descr") . "\n";
        
        foreach my $entry (values %query_to_top_hit) {
            
            my $db_id = $entry->{db_id};
            my $db_match_len = $entry->{db_match_len};
            my $db_seq_len = $seq_lengths{$db_id} or die "Error, no length found for $db_id, with hit: " . Dumper($entry);
            
            my $percent_length_matched = sprintf("%.2f", $db_match_len / $db_seq_len * 100);
            
            my $line = $entry->{line};
            my $header = $seq_headers{$db_id};
            
            print $ofh join("\t", $line, $db_match_len, $percent_length_matched, $header) . "\n";
            
            $entry->{db_hit_pct_cov} = $percent_length_matched;


            if ( ! exists $db_id_to_greatest_pct_cov{$db_id} ) {
                $db_id_to_greatest_pct_cov{$db_id} = $entry;
            }
            else {
                my $prev_entry = $db_id_to_greatest_pct_cov{$db_id};
                if ($percent_length_matched > $prev_entry->{db_hit_pct_cov}
                    ||
                    ($percent_length_matched == $prev_entry->{db_hit_pct_cov}
                     &&
                     $entry->{bitscore} > $prev_entry->{bitscore}) 
                    ) {
                    $db_id_to_greatest_pct_cov{$db_id} = $entry;
                }
            }

            $counter++;
            if ($counter % 100 == 0) {
                print STDERR "\r[$counter]    " if $verbose;
            }

        }
        close $ofh;
        
    }

    $counter = 0;
    print STDERR "\n" if $verbose;
    
    ## histogram summary

    my @bins = qw(10 20 30 40 50 60 70 80 90 100);
    my %bin_counts;

    open (my $ofh, ">$output_prefix.hist") or die "Error, cannot write to $output_prefix.hist";
    open (my $list_ofh, ">$output_prefix.hist.list") or die $!;
    {
                

        foreach my $entry (values %db_id_to_greatest_pct_cov) {

            my $pct_cov = $entry->{db_hit_pct_cov};
            
            my $prev_bin = 0;
            foreach my $bin (@bins) {
                if ($pct_cov > $prev_bin && $pct_cov <= $bin) {
                    $bin_counts{$bin}++;
                    print $list_ofh join("\t", "Bin_$bin", $entry->{line}) . "\n";
                }
                $prev_bin = $bin;
            }
            

        }
    }
    close $list_ofh;

    ## Report counts per bin
    print "#hit_pct_cov_bin\tcount_in_bin\t>bin_below\n";
    print $ofh "#hit_pct_cov_bin\tcount_in_bin\t>bin_below\n";

    my $cumul = 0;
    foreach my $bin (reverse(@bins)) {
        my $count = $bin_counts{$bin} || 0;
        $cumul += $count;
        print join("\t", $bin, $count, $cumul) . "\n";
        print $ofh join("\t", $bin, $count, $cumul) . "\n";
        
    }
    close $ofh;
    
    
    exit(0);
    

}
