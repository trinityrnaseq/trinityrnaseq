#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use SAM_reader;
use SAM_entry;
use Fasta_reader;

my $usage = "usage: $0 coord_sorted.sam genome.fa [SS_lib_type]\n\n";

my $sam_file = $ARGV[0] or die $usage;
my $genome_fa = $ARGV[1] or die $usage;
my $SS_lib_type = $ARGV[2];

main: {

	my $sam_reader = new SAM_reader($sam_file);

    my $fasta_reader = new Fasta_reader($genome_fa);
    my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();


	my $genome_acc = "";
	my $genome_seq = "";


	
	while ($sam_reader->has_next()) {
		
		my $num_introns = 0;
		
		my $sam_entry = $sam_reader->get_next();

        my $read_name = $sam_entry->reconstruct_full_read_name();
        
		my $cigar = $sam_entry->get_cigar_alignment();
		my $scaffold = $sam_entry->get_scaffold_name();

        if ($scaffold eq "*") { 
            # unaligned read
            next;
        }


		#print "$scaffold\t$cigar\n";
		
		unless ($genome_acc eq $scaffold) {
			$genome_acc = $scaffold;
			$genome_seq = $genome_seq = $genome_seqs{$genome_acc} or die "Error, no genome seq for acc: $genome_acc";
		}
        
		my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();
		
		my $strand = $sam_entry->get_query_strand();

		my $transcribed_strand = "?";
        if ($SS_lib_type) {
            $transcribed_strand = $sam_entry->get_query_transcribed_strand($SS_lib_type);
        }
		
		my $report_txt = "";

		if (scalar @$genome_coords_aref > 1) {

			my @genome_coords = @$genome_coords_aref;
			## get the introns:
			
			for (my $i = 1; $i <= $#genome_coords; $i++) {
				my $prev_coordset = $genome_coords[$i-1];	
				my $curr_coordset = $genome_coords[$i];

				my ($a_lend, $a_rend) = @$prev_coordset;
				my ($b_lend, $b_rend) = @$curr_coordset;
				
                my $intron_lend = $a_rend + 1;
                my $intron_rend = $b_lend - 1;
				
				my $intron_length = $b_lend - $a_rend - 1;
				if ($intron_length < 0) {
					die "Error, intron length invalid";
				}
				
				if ($intron_length >= 20) {
					my $intron_seq = substr($genome_seq, $a_rend +1 -1, $intron_length);
					
					my $left_intron = substr($intron_seq, 0, 2);
					my $right_intron = substr($intron_seq, -2);
					$num_introns++;
					
                    my $intron_dinucs = "$left_intron..$right_intron";
                    if ($transcribed_strand eq "?") {

                        ## make an educated guess
                        if ($intron_dinucs eq "GT..AG"
                            ||
                            $intron_dinucs eq "GC..AG"
                            ||
                            $intron_dinucs eq "AT..AC") {
                            $transcribed_strand = "+";
                        }
                        elsif ($intron_dinucs eq "CT..AC"
                               ||
                               $intron_dinucs eq "CT..GC"
                               || 
                               $intron_dinucs eq "GT..AT") {
                            $transcribed_strand =  '-';
                        }
                    }

                    
                    $report_txt .= join("\t", $read_name,
                                        "$scaffold",
                                        "$intron_lend-$intron_rend",
                                        "$transcribed_strand",
                                        "$intron_dinucs") . "\n";
				}
			}
			if ($num_introns >= 1) {
				print "$report_txt\n"; # spacer between SAM entries
			}
		}
	}


	exit(0);
}
