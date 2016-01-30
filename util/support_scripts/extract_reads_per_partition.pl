#!/usr/bin/env perl

use strict;
use warnings;
use File::Path;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use File::Basename;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use Nuc_translator;
use SAM_reader;
use SAM_entry;


my $usage = <<__EOUSAGE__;

################################################################
#
# Required:
#
# --partitions_gff <string>     list of partitions in gff format.
# --coord_sorted_SAM <string>   coordinate-sorted SAM file
#
# Options:
#
# --SS_lib_type <string>        [SS_lib_type=F,R,FR,RF]
#
# --parts_per_directory <int>  default: 100
# --min_reads_per_partition <int>   default: 10
#
#################################################################

__EOUSAGE__

    ;


my $partitions;
my $alignments_sam;
my $SS_lib_type;

my $PARTS_PER_DIR = 100;

my $MIN_READS_PER_PARTITION = 10;


&GetOptions (
             'partitions_gff=s' => \$partitions,
             'coord_sorted_SAM=s' => \$alignments_sam,
             'SS_lib_type=s' => \$SS_lib_type,
             'parts_per_directory=i' => \$PARTS_PER_DIR,
             'min_reads_per_partition=i' => \$MIN_READS_PER_PARTITION,
             );

unless ($partitions && $alignments_sam) {
    die $usage;
}


main: {


	my $partitions_dir = "Dir_". basename($partitions);
	unless (-d $partitions_dir) {
        mkdir ($partitions_dir) or die "Error, cannot mkdir $partitions_dir";
    }
	open (my $track_fh, ">$partitions_dir.listing") or die $!;
	
	my %scaff_to_partitions = &parse_partitions($partitions);

	my @ordered_partitions;
	my $current_partition = undef;
	
	my $ofh;

	my $sam_ofh;
	
	my $partition_counter = 0;

    my $part_file = "";
    my $sam_part_file = "";
    my $read_counter = 0;

    my $current_scaff = "";

    my $sam_reader = new SAM_reader($alignments_sam);
    while (my $sam_entry = $sam_reader->get_next()) {
        
		my $acc = $sam_entry->reconstruct_full_read_name();
		my $scaff = $sam_entry->get_scaffold_name();
                next if $scaff eq '*';
		my $seq = $sam_entry->get_sequence();
		my $start = $sam_entry->get_aligned_position();
        
        my $read_name = $sam_entry->get_read_name(); # raw from sam file
        if ($acc !~ /\/[12]$/ && $read_name =~ /\/[12]$/) {
            $acc = $read_name;
        }
        
        
		if (! exists $scaff_to_partitions{$scaff}) { 
			# no partitions... should explore why this is.
			print STDERR "-warning, no read partitions defined for scaffold: $scaff\n";
			next; 
		} 
		
				
		my $aligned_strand = $sam_entry->get_query_strand();
		my $opposite_strand = ($aligned_strand eq '+') ? '-' : '+';
		
		if ($aligned_strand eq '-') {
			# restore to actual sequenced bases
			$seq = &reverse_complement($seq);
		}
		
		if ($SS_lib_type) {
			## got SS data
			
			my $transcribed_orient;

            if (! $sam_entry->is_paired()) {
				if ($SS_lib_type !~ /^(F|R)$/) {
					confess "Error, read is not paired but SS_lib_type set to paired: $SS_lib_type\nread:\n$_";
				}
				
				if ($SS_lib_type eq "R") {
					$seq = &reverse_complement($seq);
				}
			}
			
			else {
				## Paired reads.
				if ($SS_lib_type !~ /^(FR|RF)$/) {
					confess "Error, read is paired but SS_lib_type set to unpaired: $SS_lib_type\nread:\n$_";
				}
                
				my $first_in_pair = $sam_entry->is_first_in_pair();
				if ( ($first_in_pair && $SS_lib_type eq "RF")
					 ||
					 ( (! $first_in_pair) && $SS_lib_type eq "FR")
					) {
					$seq = &reverse_complement($seq);
				}
			}
		}
		
	
		my $new_partition_flag = 0;

		## prime ordered partitions if first entry or if switching scaffolds.
        if ($scaff ne $current_scaff) {
                        
            $current_scaff = $scaff;
            
            @ordered_partitions = @{$scaff_to_partitions{$scaff}};                
            $current_partition = shift @ordered_partitions;
            $partition_counter = 0;

            $new_partition_flag = 1;
            

        }
        elsif ($current_partition && $start > $current_partition->{rend}) {
            $current_partition = shift @ordered_partitions;
            $partition_counter++;

            $new_partition_flag = 1;

        }
        
        if ($new_partition_flag) {
        
            close $ofh if $ofh;
			$ofh = undef;
            
			close $sam_ofh if $sam_ofh;
			$sam_ofh = undef;
			
            if ($read_counter < $MIN_READS_PER_PARTITION) {
                # delete these read files.
                #print STDERR "-- too few reads ($read_counter), removing partition: $part_file\n";
                unlink($part_file, $sam_part_file);
            }
            
            $read_counter = 0;
            
        }
        
        
        ## check to see if we're in a partition
		if (defined($current_partition) && $start >= $current_partition->{lend} && $start <= $current_partition->{rend}) {
			# may need to start a new ofh for this partition if not already established.
			unless ($ofh) {
				# create new one.
				my $file_part_count = int($partition_counter/$PARTS_PER_DIR);
				my $outdir = "$partitions_dir/" . $current_partition->{scaff} . "/$file_part_count";
				$outdir =~ s/[\;\|]/_/g;
                
				mkpath($outdir) if (! -d $outdir);
				unless (-d $outdir) {
					die "Error, cannot mkdpath $outdir";
				}
				
                $part_file = "$outdir/" . join("_", $current_partition->{lend}, $current_partition->{rend}) . ".trinity.reads";
				open ($ofh, ">$part_file") or die "Error, cannot write ot $part_file";
				print STDERR "-writing to $part_file\n";
				
				print $track_fh join("\t", $scaff, $current_partition->{lend}, $current_partition->{rend}, $part_file) . "\n";
				
                $sam_part_file = "$outdir/" . join("_", $current_partition->{lend}, $current_partition->{rend}) . ".sam";
				open ($sam_ofh, ">$sam_part_file") or die "Error, cannot open $sam_part_file";
				

			}
			# write to partition
			print $ofh ">$acc\n$seq\n";
			print $sam_ofh join("\t", $sam_entry->get_fields()) . "\n";;
            $read_counter++;

		}
        
	}
	close $track_fh;
	
	close $ofh if $ofh;
	close $sam_ofh if $sam_ofh;

	exit(0);
}



####
sub parse_partitions {
	my ($partitions_file) = @_;
	
	my %scaff_to_parts;

	print STDERR "// parsing paritions.\n";
	my $counter = 0;

	open (my $fh, $partitions_file) or die "Error, cannot open file $partitions_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { next; }
		unless (/\w/) { next; }

		$counter++;
		print STDERR "\r[$counter]  " if $counter % 100 == 0;
		
		my @x = split(/\t/);

		my $scaff = $x[0];
		my $lend = $x[3];
		my $rend = $x[4];
		my $orient = $x[6];
		
		push (@{$scaff_to_parts{$scaff}}, { scaff => $scaff,
											lend => $lend,
											rend => $rend, } );
		
	}
	print STDERR "\r[$counter]  ";
	
	close $fh;
	
	# should be sorted, but let's just be sure:
	foreach my $scaff (keys %scaff_to_parts) {
		@{$scaff_to_parts{$scaff}} = sort {$a->{lend}<=>$b->{lend}} @{$scaff_to_parts{$scaff}};
	}
	
	return(%scaff_to_parts);
}
			
	
