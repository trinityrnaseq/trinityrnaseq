#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use SAM_entry;

my $help_flag;

my $usage = <<__EOUSAGE__;

###########################################
#
#  Required:
#
# --bam <string>     bam file
#
# Optional:
#
# --max_cov <int>    default: 200
#
###########################################

__EOUSAGE__

    ;

my $bam_file;
my $max_cov = 200;

&GetOptions ( 'h' => \$help_flag,
              'bam=s' => \$bam_file,
              'max_cov=i' => \$max_cov
    );


if ($help_flag) {
    die $usage;
}

unless ($bam_file) {
    die $usage;
}

my $MAX_COV_DIST = 1e6;

main: {

    my $curr_coord = 0;
    my @cov_array;
    my %retain;
    my $prev_scaffold_name = "";
    
    open(my $fh, "samtools view -h $bam_file | ") or die "Error, cannot view bam file: $bam_file";
    while(<$fh>) {
        if (/^\@/) {
            # header line
            print $_;
            next;
        }
        my $sam_line = $_;
        my $sam_entry = new SAM_entry($sam_line);
        my $scaffold_name = $sam_entry->get_scaffold_name();
        if ($scaffold_name ne $prev_scaffold_name) {
            print STDERR "-processing $scaffold_name\n";
            ## reinit
            $curr_coord = 0;
            @cov_array = ();
            %retain = ();
        }

        my $mate_token;
        my $mate_scaffold_pos;
        if ($sam_entry->is_paired() &&
            $sam_entry->is_proper_pair() &&
            (! $sam_entry->is_mate_unmapped()) ) {

            $mate_scaffold_pos = $sam_entry->get_mate_scaffold_position();
            my $read_name = $sam_entry->get_core_read_name();
            $mate_token = join("$;", $read_name, $mate_scaffold_pos);

        }
        
        
        my $aligned_pos = $sam_entry->get_aligned_position();
        if ($curr_coord < $aligned_pos) {
            # shift it all over.
            my $delta = $aligned_pos - $curr_coord;
            if ($#cov_array > $delta) {
                @cov_array = @cov_array[$delta .. $#cov_array];
            }
            else {
                @cov_array = ();
            }
            $curr_coord = $aligned_pos;
        }
        # fill the cov array:
        my @align_coords = $sam_entry->get_alignment_coords();
        my $genome_align_coords = $align_coords[0];
        foreach my $coordset (@$genome_align_coords) {
            &add_coverage(\@cov_array, $curr_coord, $coordset);
        }
        
        my $pt_cov = $cov_array[0];
        
        #print STDERR "$curr_coord\t$pt_cov\n";
        
        if ($mate_token) {
            if ($retain{$mate_token}) {
                # force accept
                delete $retain{$mate_token};
                print $sam_line;
                next;
            }
            else {
                # was the mate already seen and not retained?
                if ($mate_scaffold_pos < $aligned_pos) {
                    next;
                }
            }
            
        }
        
        if ($pt_cov <= $max_cov) {
            # selected.
            print $sam_line;
            if ($mate_token) {
                # be sure to retain the mate too:
                $retain{$mate_token} = 1;
            }
        }
                
        $prev_scaffold_name = $scaffold_name;
        
    }
}

exit(0);

####
sub add_coverage {
    my ($cov_array_aref, $start_coord, $coordset) = @_;

    my ($lend, $rend) = @$coordset;
    $lend -= $start_coord;
    $rend -= $start_coord;

    if ($rend > $MAX_COV_DIST) {
        print STDERR "-exceeding max cov dist: start = $start_coord, lend=$lend, rend=$rend.  Skipping...\n";
        return;
    }

    for (my $i = $lend; $i <= $rend; $i++) {
        $cov_array_aref->[$i]++;
    }
}

