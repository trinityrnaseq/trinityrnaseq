#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");

use PSL_parser;
require "overlapping_nucs.ph";


# note, align to genes

my $usage = "usage: $0 blat_output.pslx [STRAND_SPECIFIC=0]\n\n";

my $blat_output = $ARGV[0] or die $usage;

my $strand_specific = $ARGV[1] || 0;


my $RESTRICTIVE_PERCENT_OVERLAP = 30; # percent of shorter length that precludes tiering
my $min_identity = 99;

my $DEBUG = 0;



main: {
    
	my %gene_to_overlapping_assemblies = &examine_blat_mappings($blat_output);

    open (my $ofh_tiers, ">$blat_output.tiers") or die $!;
    open (my $ofh_hist, ">$blat_output.tiers.hist") or die $!;
    
    
    my %tier_hist;
	
	foreach my $gene (keys %gene_to_overlapping_assemblies) {
		
		print "// analyzing gene: $gene\n" if $DEBUG;
		
		my @matches = @{$gene_to_overlapping_assemblies{$gene}};
		
		@matches = sort {$a->{score}<=>$b->{score}} @matches;


		my @tiers;
		my $match = pop @matches;
		push (@tiers, [$match]);
		
		while (@matches) {
			my $match = pop @matches;
			my $placed_flag = 0;

			print "\t-analyzing match: " . join("\t", $match->{acc}, $match->{gene_lend}, $match->{gene_rend}) . "\n" if $DEBUG;
			
			

			foreach my $tier (@tiers) {
				
				if (! &restrictive_overlap($match, $tier)) {
					push (@$tier, $match);
					$placed_flag = 1;
					print "* match placed\n" if $DEBUG;
					last;
				}
			}
			unless ($placed_flag) {
				print "* match NOT placed\n" if $DEBUG;
                print "** adding match to new tier.\n" if $DEBUG;
                push (@tiers, [$match]);
            }
		}
		
		## report matches that were tiered.
		
		my $tier_counter = 0;
		foreach my $tier_aref (@tiers) {
			
			$tier_counter++;
			
			my @tier = @$tier_aref;
			@tier = sort {$a->{gene_lend}<=>$b->{gene_lend}} @tier;
			
			my $gene_token = $gene;
			# $gene_token =~ s/,anti//;
			
			print $ofh_tiers "Tier[$tier_counter]\t$gene";
			
            $tier_hist{$tier_counter}++;
            
            foreach my $match (@tier) {
				my $acc = $match->{acc};
				my $gene_lend = $match->{gene_lend};
				my $gene_rend = $match->{gene_rend};
				my $per_id = $match->{per_id};
				my $trans_end5 = $match->{trans_end5};
				my $trans_end3 = $match->{trans_end3};
				my $total_gaps = $match->{total_gaps};
							
				my $gene_len = $match->{gene_len};
				my $trans_len = $match->{trans_len};
				
				my $score = $match->{score};

				my $matches = $match->{matches};
				my $mismatches = $match->{mismatches};
				
				print $ofh_tiers "\t$acc\{S[$score],G:$gene_lend-$gene_rend:$gene_len,T:$trans_end5-$trans_end3:$trans_len,$per_id,M:$matches,N:$mismatches,gaps:$total_gaps}";
			}
			print $ofh_tiers "\n";
			            			
		}
	}
	
    close $ofh_tiers;

    ## write histogram
    my @tiers = sort {$a<=>$b} keys %tier_hist;
    foreach my $tier (@tiers) {
        my $count = $tier_hist{$tier};
        print $ofh_hist join("\t", $tier, $count) . "\n";
    }
    close $ofh_hist;
    
	exit(0);
	
	
}


####
sub examine_blat_mappings {
	my ($blat_output) = @_;
	
	my %mappings;
	
	my $psl_reader = new PSL_parser($blat_output);
	
	while (my $psl_entry = $psl_reader->get_next()) {
		
		my $trans_assembly = $psl_entry->get_Q_name();
		my $gene_id = $psl_entry->get_T_name();
		my $per_id = $psl_entry->get_per_id();
		
		my $trans_len = $psl_entry->get_Q_size();
		my $gene_len = $psl_entry->get_T_size();
		

		my ($trans_end5, $trans_end3) = $psl_entry->get_Q_span();
		my ($gene_end5, $gene_end3) = $psl_entry->get_T_span();
		
		my $strand = $psl_entry->get_strand();
		my $total_gaps = $psl_entry->get_T_gap_bases() + $psl_entry->get_Q_gap_bases();
		my $num_matches = $psl_entry->get_match_count();
		my $num_mismatches = $psl_entry->get_mismatch_count();
		

		my $percent_gaps = ($total_gaps) / ($num_matches + $num_mismatches) * 100; 
		
		if ($per_id < $min_identity) { next; }

				
		my ($gene_lend, $gene_rend) = sort {$a<=>$b} ($gene_end5, $gene_end3);
	
		my $score = (5 * $psl_entry->get_match_count()) - (4 * $psl_entry->get_mismatch_count()) - (1 * $total_gaps);
		
		my $token = ($strand eq '-' && $strand_specific) ? ",anti" : "";
		
		push (@{$mappings{"$gene_id$token"}}, { acc => "$trans_assembly", 
												
												gene_lend => $gene_lend,
												gene_rend => $gene_rend,
												gene_len => $gene_len,
												
												per_id => $per_id,
												score => $score,
												
												trans_end5 => $trans_end5,
												trans_end3 => $trans_end3,
												trans_len => $trans_len,
												
												matches => $psl_entry->get_match_count(),
												mismatches => $psl_entry->get_mismatch_count(),
												total_gaps => $total_gaps,

											} );
		
		
	}
	
	if ($DEBUG) {
		foreach my $gene (keys %mappings) {
			print "gene: $gene\n";
			foreach my $mapping (@{$mappings{$gene}}) {
				print "\t" . join ("\t", $mapping->{acc}, $mapping->{gene_lend}, $mapping->{gene_rend}) . "\n";
			}
		}
	}
	
	return(%mappings);
}






####
sub restrictive_overlap {
	my ($match, $tier_aref) = @_;

	my ($match_lend, $match_rend) = ($match->{gene_lend}, $match->{gene_rend});

	my $match_len = $match_rend - $match_lend + 1;

	foreach my $ele (@$tier_aref) {
		
		my ($ele_lend, $ele_rend) = ($ele->{gene_lend}, $ele->{gene_rend});
		
		## check for encapsulations:
		if (&coordset_A_encapsulates_B( [$match_lend, $match_rend], [$ele_lend, $ele_rend] ) 
			||
			&coordset_A_encapsulates_B( [$ele_lend, $ele_rend], [$match_lend, $match_rend] ) ) {
			return(1);
		}

		my $ele_len = $ele_rend - $ele_lend + 1;
		
		if (&coordsets_overlap ([$ele_lend, $ele_rend], [$match_lend, $match_rend] ) ) {
		
			my $overlap_len = &nucs_in_common($match_lend, $match_rend, $ele_lend, $ele_rend);
			
			if ( $overlap_len / $match_len * 100 >= $RESTRICTIVE_PERCENT_OVERLAP
				 ||
				 $overlap_len / $ele_len * 100 >= $RESTRICTIVE_PERCENT_OVERLAP) {
				return(1);
			}
		}
	}
	
	return(0); # no restrictive overlap
}


		
