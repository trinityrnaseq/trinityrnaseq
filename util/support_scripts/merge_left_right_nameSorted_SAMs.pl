#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SAM_reader;
use SAM_entry;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

##########################################################
#
#  --left_sam    left sam file, read name and contig name -sorted (sort -k1,1 -k3,3 file.sam > file.sorted.sam)
#  --right_sam   right sam file, read name and contig name -sorted
#  
#  -D            maximum distance between pairs (default: 2000)
#  
#  
#
###########################################################

_EOUSAGE_

        ;


my ($left_sam_file, $right_sam_file);
my $max_dist_between_pairs = 2000;


&GetOptions( 'left_sam=s' => \$left_sam_file,
			 'right_sam=s' => \$right_sam_file,
			 'D=i' => \$max_dist_between_pairs,
	);

unless ($left_sam_file && $right_sam_file) {
	die $usage;
}



main: {

	my $left_sam_reader = new SAM_reader($left_sam_file);
	
	my $right_sam_reader = new SAM_reader($right_sam_file);


	my $iter_counter = 0;
	
	while ($left_sam_reader->has_next() && $right_sam_reader->has_next()) {

		$iter_counter++;
		print STDERR "\r[$iter_counter]  " if ($iter_counter % 10000 == 0);
		
		## process Left read, unpaired
		if (
            $left_sam_reader->preview_next()->get_core_read_name() lt $right_sam_reader->preview_next()->get_core_read_name()
                        
            ) {
			
			my @left_entries;
			
			my $next_sam = $left_sam_reader->get_next();
			push (@left_entries, $next_sam);
			
			while ($left_sam_reader->has_next() 
				   && 
				   $left_sam_reader->preview_next()->get_core_read_name() eq $left_entries[$#left_entries]->get_core_read_name()
                   ) {
				
				$next_sam = $left_sam_reader->get_next();
				push (@left_entries, $next_sam);
			}
            
            while (@left_entries) {
                
                $next_sam = shift @left_entries; # take the one with the highest mapping quality
                
                $next_sam->set_first_in_pair(1);
                $next_sam->set_mate_unmapped(1);
                $next_sam->set_paired(1);
                
                print $next_sam->toString(1) . "\n";
            }
        }
		
		## process right read, unpaired
		elsif (
               $left_sam_reader->preview_next()->get_core_read_name() gt $right_sam_reader->preview_next()->get_core_read_name()
               
            ) {
			
            
			my @right_entries;
			
			my $next_sam = $right_sam_reader->get_next();
			push (@right_entries, $next_sam);
			
			while ($right_sam_reader->has_next()
				   &&
				   $right_sam_reader->preview_next()->get_core_read_name() eq $right_entries[$#right_entries]->get_core_read_name()
                                      
                ) {
				
				$next_sam = $right_sam_reader->get_next();
				push (@right_entries, $next_sam);
			}
			
			
            while (@right_entries) {
			
                $next_sam = shift @right_entries;
			
                $next_sam->set_second_in_pair(1);
                $next_sam->set_mate_unmapped(1);
                $next_sam->set_paired(1);
                
                print $next_sam->toString(1) . "\n";
            }
        }
		
		## process paired-aligned reads
		elsif (
               $left_sam_reader->preview_next()->get_core_read_name() eq $right_sam_reader->preview_next()->get_core_read_name()
             
            ) {  # must be the same
			
			## process left entries
			my @left_entries = $left_sam_reader->get_next();
			while ($left_sam_reader->has_next() 
				   && 
				   $left_sam_reader->preview_next()->get_core_read_name() eq $left_entries[$#left_entries]->get_core_read_name()
               
                ) {
				
				push (@left_entries, $left_sam_reader->get_next());
			}
			
			## process right entries
			my @right_entries = $right_sam_reader->get_next();
			while ($right_sam_reader->has_next()
				   &&
				   $right_sam_reader->preview_next()->get_core_read_name()  eq $right_entries[$#right_entries]->get_core_read_name()
                   ) {
				
                push (@right_entries, $right_sam_reader->get_next());
			}
			
			
			&identify_mapped_pairs(\@left_entries, \@right_entries, $max_dist_between_pairs);
			
			
		}

		elsif ($left_sam_reader->has_next() && $right_sam_reader->has_next()) {
			
            ## Shouldn't ever get here. Fatal error.

            my $left_read_name = $left_sam_reader->preview_next()->get_core_read_name();
			my $right_read_name = $right_sam_reader->preview_next()->get_core_read_name();

			confess("Error, left and right read names do not appear comparable: left=$left_read_name, right=$right_read_name");
		}
		
	}
	

	## Finish printing the last ones.
	
	print STDERR "\nreporting remaining left reads\n";
	
	while ($left_sam_reader->has_next()) {

        my $next_sam = $left_sam_reader->get_next();
        

        $next_sam->set_first_in_pair(1);
        $next_sam->set_mate_unmapped(1);
        $next_sam->set_paired(1);
        
        print $next_sam->toString(1) . "\n";
    }
    	
	print STDERR "\nreporting remaining right reads\n";
	while ($right_sam_reader->has_next()) {

        my $next_sam = $right_sam_reader->get_next();
        $next_sam->set_second_in_pair(1);
        $next_sam->set_mate_unmapped(1);
        $next_sam->set_paired(1);
        
        print $next_sam->toString(1) . "\n";
    }
    	
	print STDERR "\nDone.\n\n";
	
	exit(0);
	
}
	

	
####
sub identify_mapped_pairs {
	my ($left_entries_aref, $right_entries_aref, $max_dist_between_pairs) = @_;

    
    ## group according to scaffold.
    
    my %scaffold_to_pairs;
    
    foreach my $left_sam_entry (@$left_entries_aref) {
        
        my $scaffold = $left_sam_entry->get_scaffold_name();
        
        push (@{$scaffold_to_pairs{$scaffold}->{left}}, $left_sam_entry);
    }

    foreach my $right_sam_entry (@$right_entries_aref) {
        
        my $scaffold = $right_sam_entry->get_scaffold_name();
        
        push (@{$scaffold_to_pairs{$scaffold}->{right}}, $right_sam_entry);
    }

    foreach my $scaffold (keys %scaffold_to_pairs) {
        
        my $left_aref = [];
        my $right_aref = [];

        if (my $aref =  $scaffold_to_pairs{$scaffold}->{left}) {
            $left_aref = $aref;
        }
        if (my $aref = $scaffold_to_pairs{$scaffold}->{right}) {
            $right_aref = $aref;
        }
        
        &identify_mapped_pairs_single_scaffold($left_aref, $right_aref, $max_dist_between_pairs);
    }
    
    return;
}

sub identify_mapped_pairs_single_scaffold {
    my ($left_entries_aref, $right_entries_aref, $max_dist_between_pairs) = @_;

	my $num_left = scalar @$left_entries_aref;
	my $num_right = scalar @$right_entries_aref;

	#print STDERR "\t$num_left, $num_right\n";
	

	# Sort by mapping quality (this puts tophat alignments ahead of blat alignments for now, since blat sam doesn't include quality values)
	@$left_entries_aref = &sort_by_score_descending(@$left_entries_aref);
	@$right_entries_aref = &sort_by_score_descending(@$right_entries_aref);
	
	
	foreach my $left_entry (@$left_entries_aref) {
		$left_entry->set_first_in_pair(1);
	}
	
	foreach my $right_entry (@$right_entries_aref) {
		$right_entry->set_second_in_pair(1);
	}

	foreach my $entry (@$left_entries_aref, @$right_entries_aref) {
		$entry->set_paired(1);
	}

	my @pairs;
	
	my @examined_left;
	my @examined_right;
	
	
	# find the pair that has the highest mapping quality
	for (my $i = 0; $i <= $#$left_entries_aref; $i++) {

		my $left_entry = $left_entries_aref->[$i];

		$examined_left[$i] = 0;

		for (my $j = 0; $j <= $#$right_entries_aref; $j++) {
			
			if (defined $examined_right[$j] && $examined_right[$j]) { 
				next;
			}
			
			$examined_right[$j] = 0;
			
			my $right_entry = $right_entries_aref->[$j];
			
			if ($left_entry->get_query_strand() ne $right_entry->get_query_strand() 
				
				&&
				
				$left_entry->get_scaffold_name() eq $right_entry->get_scaffold_name()
				
				&&
				
				abs ($left_entry->get_aligned_position() - $right_entry->get_aligned_position()) <= $max_dist_between_pairs
				
			) {
				
				# candidate proper mapped pair.  double check orientation.
				
				if ( ($left_entry->get_query_strand eq '+'
					  &&
					  $left_entry->get_aligned_position() < $right_entry->get_aligned_position()
					  )
					 
					 ||
					 
					 ($left_entry->get_query_strand eq '-'
					  && 
					  $left_entry->get_aligned_position() > $right_entry->get_aligned_position() )
					) 
				{
					
					## excellent proper mapped pairs
					push (@pairs, [$left_entry, $right_entry]);
					$examined_left[$i] = 1;
					$examined_right[$j] = 1;
					last;
				}
			}
		}
				
	}


	my $num_reads_reported = 0;

	if (@pairs) {
		
		foreach my $pair (@pairs) {
			my ($left_entry, $right_entry)  = @$pair;

			$left_entry->set_mate_strand($right_entry->get_query_strand() );
			$right_entry->set_mate_strand($left_entry->get_query_strand() );

			$left_entry->set_mate_scaffold_name('=');
			$right_entry->set_mate_scaffold_name('=');
			
			$left_entry->set_mate_scaffold_position( $right_entry->get_aligned_position() );
			$right_entry->set_mate_scaffold_position( $left_entry->get_aligned_position() );
			
			$left_entry->set_proper_pair(1);
			$right_entry->set_proper_pair(1);
			
            my $genome_span = &get_genome_span_from_pair($left_entry, $right_entry);
            
            # set template length to span on the genome (reasonable option, given rna-seq data)
            $left_entry->{_fields}->[8] = $genome_span;  ## hacky way of doing it
            $right_entry->{_fields}->[8] = $genome_span;
            
			print $left_entry->toString(1) . "\n"
				. $right_entry->toString(1) . "\n";
			
			$num_reads_reported++;
			
		}
	}

	## report each of the remaining alignments:
			
    for (my $i = 0; $i <= $#$left_entries_aref; $i++) {
        if (! $examined_left[$i]) {
            my $entry = $left_entries_aref->[$i];
            $entry->set_mate_unmapped(1);
            print $entry->toString(1) . "\n";
            
                        
        }
    }
			
	
    for (my $j = 0; $j <= $#$right_entries_aref; $j++) {
        if (! $examined_right[$j]) {
            my $entry = $right_entries_aref->[$j];
            $entry->set_mate_unmapped(1);
            print $entry->toString(1) . "\n";
            
        }
    }
		
	return;
	
}


####
sub sort_by_score_descending {
	my @entries = @_;

	my @structs;

	foreach my $entry (@entries) {
		
		my $score = &get_score($entry);
		
		push (@structs, { entry => $entry,
						  score => $score,
			  } );
		
	}

	@structs = reverse sort {$a->{score}<=>$b->{score}} @structs;

	my @ret_entries;
	foreach my $struct (@structs) {
		push (@ret_entries, $struct->{entry});
	}

	return(@ret_entries);
}

####
sub get_score {
	my ($entry) = @_;

	my $score = $entry->get_mapping_quality();
		
	## use alignment score (blat) if available.
	my $text = $entry->toString();
	if ($text =~ /AS:i:(\d+)/) {
		my $align_score = $1;
		$score = $align_score;
	}
	
	return($score);
}

####
sub sort_pairs_by_score_descending {
	my (@pairs) = @_;

	my @structs;

	foreach my $pair (@pairs) {
		
		my ($left, $right) = @$pair;

		my $score = &get_score($left) + &get_score($right);
		
		push (@structs, { pair => $pair,
						  score => $score,
			  });
	}

	@structs = reverse sort { $a->{score}<=>$b->{score}} @structs;

	my @ret_pairs;

	foreach my $struct (@structs) {
		my $pair = $struct->{pair};
		push (@ret_pairs, $pair);
	}


	return(@ret_pairs);
}


####
sub get_genome_span_from_pair {
    my ($left_entry, $right_entry) = @_;

    my @coords = ($left_entry->get_genome_span(), $right_entry->get_genome_span());
    @coords = sort {$a<=>$b} @coords;

    my $lend = shift @coords;
    my $rend = pop @coords;

    my $length = $rend - $lend + 1;

    return($length);
}
