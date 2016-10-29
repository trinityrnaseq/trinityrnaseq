#!/usr/bin/env perl

# NOTE FOR WEB GUIDE http://trinityrnaseq.sourceforge.net/genome_guided_trinity.html
# % samtools view gsnap_out/gsnap.coordSorted.bam > gsnap.coordSorted.sam
# change to run samtools view with -F4 or -F12 if you are plannign to use --no_single

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use SAM_reader;
use SAM_entry;
use COMMON;

use Getopt::Long qw(:config no_ignore_case bundling);

use Carp;
use Data::Dumper;
use Cwd;

$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script


my $usage = <<_EOUSAGE_;

#######################################################################
#
# Required:
#
# --sam <string>                SAM-formatted file
#
# Optional:
#
# --no_single                           no single reads
#
# --min_insert_size <int>               minimum distance between proper pair's starting positions. (default: 100)
# --max_insert_size <int>               maximum distance between proper pair's starting positions. (default: 500)
#
# --sort_buffer <string>                default '10G', amount of RAM to allocate to sorting
#
# -d                             debug mode
#
#######################################################################


_EOUSAGE_

    ;

my $sam_file;
my $full_flag = 0;
my $MAX_INSERT_SIZE = 500;
my $MIN_INSERT_SIZE = 100;
my $DEBUG = 0;
my $help_flag = 0;
my $NO_SINGLE;
my $CPU = 1;
my $sort_buffer = '10G';
&GetOptions( 
             'h' => \$help_flag,

             'sam=s' => \$sam_file,
             'CPU=i' => \$CPU,
             'sort_buffer=s' => \$sort_buffer,
             'no_single' => \$NO_SINGLE,
             'max_insert_size=i' => \$MAX_INSERT_SIZE,
             'min_insert_size=i' => \$MIN_INSERT_SIZE,
             
             'd' => \$DEBUG,

             );


if ($help_flag) {
    die $usage;
}


unless ($sam_file) {
    die $usage;
}

my $sort_exec = &COMMON::get_sort_exec($CPU);


main: {
	

	my @frags;

    my $read_coords_file = "$sam_file.read_coords";
    my $read_coords_checkpoint = "$read_coords_file.ok";
    
    &extract_read_coords($read_coords_file) unless (-s $read_coords_checkpoint);
    &process_cmd("touch $read_coords_checkpoint");
    
    my $pair_coords_file = "$sam_file.frag_coords";
    &extract_frag_coords($read_coords_file, $pair_coords_file);
            
    
    exit(0);



}

####
sub extract_read_coords {
    my ($read_coords_file) = @_;
    # AP: I think this is very expensive (more than 20' for 12G SAM)
    # option 1) the same information could be derived by parsing the BAM file
    # and cut -f 1,2,3,7,8,4,9,10 -> there will be two identical (1) fields, if paired the one with (9)>0 is the /1 one||the one on the plus strand (from (2)) could be assigned as the /1 one
    # generally, i don't see the need to have to create a SAM file in the first place... 
    # another benefit of samtools view bam <region> is that we could extract each scaffold separately (would help with sort...)
    # option 2) use samtools depth instead of wig file???
    
    ## Every read processed.
    
    print STDERR "-extracting read coordinates from $sam_file into $read_coords_file\n\n";
    
    my $sam_reader = new SAM_reader($sam_file);
    
    open (my $ofh, ">$read_coords_file") or die "Error, cannot write to file $read_coords_file";
    
    while ($sam_reader->has_next()) {
        
        my $sam_entry = $sam_reader->get_next();
        
        # unless ($sam_entry->get_mate_scaffold_name() eq "=" || $sam_entry->get_mate_scaffold_name() eq $sam_entry->get_scaffold_name()) { next; }
        # commenting out above - just describe the reads, let the other routine handle the fragment definition.


        eval {
            my $scaffold = $sam_entry->get_scaffold_name();
            my $core_read_name = $sam_entry->get_core_read_name();
            my $read_name = $sam_entry->get_read_name();
            my $full_read_name = $sam_entry->reconstruct_full_read_name();
            
            #print "read_name: $read_name, full_read_name: $full_read_name\n";
            
            
            my $pair_side = ".";
            if ($full_read_name =~ m|/([12])$|) {
                $pair_side = $1;
            }
            elsif ($read_name =~ m|/([12])$|) {
                $pair_side = $1;
            }
            
            my $mate_scaff_pos = $sam_entry->get_mate_scaffold_position();
            my ($read_start, $read_end) = $sam_entry->get_genome_span();
            
            
            print $ofh join("\t", $scaffold, $core_read_name, $pair_side, $read_start, $read_end) . "\n" if $scaffold ne '*' && $read_start && $read_end;
        };

        if ($@) {
            print STDERR "******\nError parsing SAM entry: " . Dumper($sam_entry) . " \n$@\n******\n\n\n";
        }
    }
    
    close $ofh;

    
    return;
        
}


####
sub extract_frag_coords {
    my ($read_coords_file, $pair_frag_coords_file) = @_;
    unless (-s "$read_coords_file.sort_by_readname" ){
        ## sort by scaffold, then by read name
        my $cmd = "$sort_exec -S$sort_buffer -T . -k1,1 -k2,2 -k4,4n $read_coords_file > $read_coords_file.sort_by_readname";
        &process_cmd($cmd);
        rename("$read_coords_file.sort_by_readname", $read_coords_file);
        # so that sort is not re-done...
        symlink(&create_full_path($read_coords_file),&create_full_path("$read_coords_file.sort_by_readname"));
    }
    
    ## define fragment pair coordinate span
    open (my $fh, "$read_coords_file") or die $!;
    
    open (my $ofh, ">$pair_frag_coords_file") or die $!;
    
    my $prev_reported_pair = "";
    my $prev_reported_single = "";
    
    my $first = <$fh>;
    chomp $first;
    while (my $second = <$fh>) {
        next if ($first =~/^\*/ && $second =~/^\*/); # both unmapped
        chomp $second;
        
        my ($scaffA, $readA, $readA_pair_side, $lendA, $rendA) = split(/\t/, $first);
        my ($scaffB, $readB, $readB_pair_side, $lendB, $rendB) = split(/\t/, $second);
        
        my $got_pair_flag = 0;
        if ($readA eq $readB 
            && $scaffA eq $scaffB 
            && $readA_pair_side ne $readB_pair_side 
            && $readA_pair_side =~ /\d/ && $readB_pair_side =~ /\d/) {
            
            my @coords = sort {$a<=>$b} ($lendA, $rendA, $lendB, $rendB);
            my $min = shift @coords;
            my $max = pop @coords;
            
            my $insert_size = $max - $min + 1;
            if ($insert_size >= $MIN_INSERT_SIZE && $insert_size <= $MAX_INSERT_SIZE && $readA ne $prev_reported_pair) {
                # treat as proper pair
                
                print $ofh join("\t", $scaffA, $readA, $min, $max) . "\n";
                $prev_reported_pair = $readA; # only one proper pair to be reported. - not random though, first one encountered.
                
            }
            else {
                # treat as unpaired reads
                unless ($NO_SINGLE) {
                    if ($prev_reported_single ne $readA) {
                        print $ofh join("\t", $scaffA, $readA . "/$readA_pair_side", $lendA, $rendA) . "\n";
                        print $ofh join("\t", $scaffB, $readB . "/$readB_pair_side", $lendB, $rendB) . "\n";
                    }
                    $prev_reported_single = $readA;
                }
            }
    
            $first = <$fh>; # prime first
            chomp $first if $first;
        }
        else {
            # not paired
            unless ($NO_SINGLE) {
                if ($prev_reported_single ne $readA) {
                    print $ofh join("\t", $scaffA, $readA . "/$readA_pair_side", $lendA, $rendA) . "\n";
                }
                $prev_reported_single = $readA; # only letting one slip through.
            }
            
            $first = $second;
            next;
        }
        
        

    }

    close $ofh;
    close $fh;


    my $cmd = "$sort_exec -S$sort_buffer -T . -k1,1 -k3,3n $pair_frag_coords_file > $pair_frag_coords_file.coord_sorted";
    &process_cmd($cmd) unless -s "$pair_frag_coords_file.coord_sorted";

    rename("$pair_frag_coords_file.coord_sorted", $pair_frag_coords_file);
        
    return;
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}

####
sub create_full_path {
    my ($path) = @_;

    unless ($path =~ /^\//) {
        $path = cwd() . "/$path";
    }

    return($path);
}
