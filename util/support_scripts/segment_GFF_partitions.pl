#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 partitions.gff clip_pts.wig\n\n";

my $partitions_gff = $ARGV[0] or die $usage;
my $clip_pts_wig = $ARGV[1] or die $usage;

my $MIN_PARTITION_SIZE = 50;

main: {

	my %scaff_to_gff = &parse_gff_entries($partitions_gff);
	my %scaff_to_clip = &parse_wig_entries($clip_pts_wig);

	
	foreach my $scaff (sort keys %scaff_to_gff) {
		
		my @gffs = @{$scaff_to_gff{$scaff}};
		@gffs = sort {$a->{lend}<=>$b->{lend}} @gffs;

		my @clips;
		if (exists $scaff_to_clip{$scaff}) {
			@clips = @{$scaff_to_clip{$scaff}};
			@clips = sort {$a<=>$b} @clips;
		}
		
		foreach my $gff (@gffs) {
			
			my ($lend, $rend) = ($gff->{lend}, $gff->{rend});
			
			my @overlapping_clips;

			if (@clips) {
				while (@clips && $clips[0] < $rend) {
					my $clip = shift @clips;
					if ($clip > $lend) {
						push (@overlapping_clips, $clip);
					}
				}
			}
			
			if (@overlapping_clips) {
				my @x = @{$gff->{line}};
				
				my $gff_line = join("\t", @x);
				
				#print "$gff_line\nHas overlapping clips: @overlapping_clips\n";
				

				my $prev_lend = $lend;
				
				foreach my $overlapping_clip (@overlapping_clips) {
					
					my $rend = $overlapping_clip - 1;
					
					my @y = @x;
					$y[3] = $prev_lend;
					$y[4] = $rend;

					my $part_length = $rend - $prev_lend + 1;
					
					if ($part_length >= $MIN_PARTITION_SIZE) {

						print join("\t", @y) . "\n";

					}
					
					$prev_lend = $rend + 1;
					
					
					#print "$gff_line\nWith clips: " . join("\t", @overlapping_clips) . "\n\n";
				}
				## report last partition of clip:
				
				$x[3] = $prev_lend;
				print join("\t", @x) . "\n";
				

				#print "\n\n";
			}

			else {
				
				## no clipping.  Report original partition.
				
				print join("\t", @{$gff->{line}}) . "\n";
			}
			

		}

	}

	exit(0);

}
			


####
sub parse_gff_entries {
	my ($gff_file) = @_;

	my %scaffold_to_gff;
	
	open (my $fh, $gff_file) or die "Error, cannot open file $gff_file";
	while (<$fh>) {
		chomp;
		
		my @x = split(/\t/);
		
		my $scaff = $x[0];
		my $lend = $x[3];
		my $rend = $x[4];
		
		my $struct = { line => [@x],
					   scaff => $scaff,
					   lend => $lend,
					   rend => $rend,
				   };

		push (@{$scaffold_to_gff{$scaff}}, $struct);
	}
	close $fh;

	return(%scaffold_to_gff);
}


#### 
sub parse_wig_entries {
	my ($jaccard_file) = @_;

	my %scaffold_to_entries;

	my $curr_scaff;
        
	open (my $fh, $jaccard_file) or die "Error, cannot open file $jaccard_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		chomp;
		if (/^variableStep chrom=(\S+)/) {
			$curr_scaff = $1;
			next;
		}
                        
		my ($coord, $val) = split(/\t/);
                
		unless (defined($coord) && defined($val)) {
                        #die "Error, line $_ lacks expected format";
			next;
		}
                
		push (@{$scaffold_to_entries{$curr_scaff}}, $coord);
	}

	close $fh;


	return(%scaffold_to_entries);
}

