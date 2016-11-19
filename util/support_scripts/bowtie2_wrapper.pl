#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case bundling);

$ENV{PATH} .= "\:$FindBin::RealBin/../trinity-plugins/rsem/sam/";  # include samtools in path, already included in rsem build.

$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script


my $output_prefix = "bowtie2";


my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  --left and --right <string>   reads
#
#        or
#
#  --single <string>             reads
#
#  Required inputs:
#   
#  --target <string>            multi-fasta file containing the target sequences (should be named {refName}.fa )
#
#  --seqType <string>           fa | fq    (fastA or fastQ format)
#
#
# Optional:
#
#  --SS_lib_type <string>       strand-specific library type:  single: F or R  paired: FR or RF
#                                   examples:  single RNA-Ligation method:  F
#                                              single dUTP method: R
#                                              paired dUTP method: RF
#
#  --num_top_hits <int>         (default: 20) 
#
#  --retain_intermediate_files     retain all the intermediate sam files produced (they take up lots of space! and there's lots of them)
#
#  --max_dist_between_pairs             default (2000) 
#
#  --just_prep_build               just prepare the bowtie-build and stop.
#
#  --output_prefix|o <string>      prefix for output filename (default: $output_prefix)
#
#  --CPU <int>                     number of threads
#
#  ## General options
#
#  Any options after '--' are passed onward to the alignments programs (except BLAT -which has certain options exposed above).
#  For example, to set the number of processors used by Bowtie to 16 then use:
#     -- -p 16
#
####################################################################################################################
#
#  Example commands:
#
#     $0 --seqType fq --left left.fq --right right.fq --target Trinity.fasta 
#
#
####################################################################################################################



_EOUSAGE_

    ;


my $help_flag;
my $target_db; 
my $left_file;
my $right_file;
my $single_file;

my $num_top_hits = 20;
my $max_dist_between_pairs = 2000;
my $seqType;
my $SS_lib_type;
my $retain_intermediate_files_flag = 0;

my $trinity_mode;
my $CPU = 1;

my $JUST_PREP_BUILD = 0;

unless (@ARGV) {
    die $usage;
}


&GetOptions ( 'h' => \$help_flag,
              
              ## required inputs
              'left=s' => \$left_file,
              'right=s' => \$right_file,
              'single=s' => \$single_file,
              
              "target=s" => \$target_db,
              "seqType=s" => \$seqType,
              
              "CPU=i" => \$CPU,
              

              ## Optional:
              "SS_lib_type=s" => \$SS_lib_type,
              
              'output_prefix|o=s' => \$output_prefix,
              
              'num_top_hits=i' => \$num_top_hits,
                            
              'max_dist_between_pairs=i' => \$max_dist_between_pairs,
              
              'retain_intermediate_files' => \$retain_intermediate_files_flag,
              
              'just_prep_build' => \$JUST_PREP_BUILD,
              
              );


if ($help_flag) { die $usage; }


unless ($target_db && -s $target_db) { 
    die $usage . "Must specify target_db and it must exist at that location";
}

unless ($JUST_PREP_BUILD || ($seqType && $seqType =~ /^(fq|fa)$/)) {
    die $usage . ", sorry do not understand seqType $seqType";
}


unless ($JUST_PREP_BUILD 
        || 
        ($left_file && $right_file)
        ||
        ($single_file)

    )  {
    die $usage . "sorry, cannot find files $left_file $right_file $single_file";
}


if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
    die "Error, SS_lib_type must be one of the following: (F, R, FR, RF)  ";
}



## check for required programs
{
    
    my @required_progs = qw(samtools bowtie2-build bowtie2);
    
    foreach my $prog (@required_progs) {
        my $path = `which $prog`;
        unless ($path =~ /^\//) {
            die "Error, path to required $prog cannot be found";
        }
    }
}


my $util_dir = "$FindBin::RealBin";

my ($start_dir, $work_dir, $num_hits);


main: {
    $start_dir = cwd();

    $left_file = &build_full_paths($left_file, $start_dir) if $left_file;
    $right_file = &build_full_paths($right_file, $start_dir) if $right_file;
    $target_db = &build_full_paths($target_db, $start_dir);
    
    $single_file = &build_full_paths($single_file, $start_dir) if $single_file;
    
    unless (-s "$target_db.fai") {
        &process_cmd("samtools faidx $target_db");
    }
    
    
    ######################################
    ## Prep the bowtie index of the target
    ######################################
    
    my $index_ext = "bt2";
    
    my @bowtie_build_files = <$target_db.*.$index_ext>;
    print Dumper(\@bowtie_build_files);
    unless (@bowtie_build_files) {
        
        print STDERR "Note - bowtie-build indices do not yet exist. Indexing genome now.\n";
        ## run bowtie-build:
        
        my $start_index_file = "$target_db._${index_ext}_idx_.start";
        if (-e $start_index_file) {
            die "Error, indexing in progress or never finished. Wait for other process to finish building the index, or delete the $start_index_file and try again";
        }
        &process_cmd("touch $start_index_file");
        # make bowtie index
        
        my $builder = "bowtie2-build";
        
        my $cmd = "$builder -q $target_db $target_db";
        &process_cmd($cmd);
        
        unlink($start_index_file);
        
    }
    
    if ($JUST_PREP_BUILD) {
        print STDERR "Just preparing build, stopping now.\n";
        exit(0);
    }
    
    my @entries;
    
    my $format = ($seqType eq "fq") ? "-q" : "-f";
    
    my $output_file = "$output_prefix.coordSorted.bam";

    my $cmd;
    if ($left_file && $right_file) {

        $cmd = "bash -c \"set -o pipefail; bowtie2 --local -k $num_top_hits --threads $CPU $format -x $max_dist_between_pairs -x $target_db -1 $left_file -2 $right_file | samtools view -F4 -Sb - | samtools sort -o - - > $output_file\" ";
        

    }
    else {
        
        $cmd = "bash -c \"set -o pipefail; bowtie2 --local -k $num_top_hits --threads $CPU $format -x $target_db -U $single_file | samtools view -F4 -Sb - | samtools sort -o - - > $output_file\" ";
        
    }
    
    &process_cmd($cmd) unless (-e "$output_file.ok");

    &process_cmd("touch $output_file.ok") unless (-e "$output_file.ok");
    
    
    my @to_delete;
    
    if ($SS_lib_type) {
        
        ## strand-specific
        ## separate the sam based on strand, and create separate bam files.  (for convenience sake)
        
        $cmd = "$util_dir/SAM_strand_separator.pl $output_file $SS_lib_type";
        &process_cmd($cmd);
        
                
        foreach my $sam_file ("$output_file.+.sam", "$output_file.-.sam") {
            
            my $bam_file = $sam_file;
            $bam_file =~ s/\.sam$/\.bam/; # add suffix below

            unless (-s $sam_file) { next; } # empty file
                
            push (@to_delete, $sam_file);
        
            $cmd = "samtools view -bt $target_db.fai $sam_file | samtools sort -o - - > $bam_file"; # .bam ext added auto
                        
            &process_cmd($cmd);
            
            $cmd = "samtools index $bam_file";
            &process_cmd($cmd);
            
        }
    }
    
        
    ## do final cleanpup of intermediate files
    unless ($retain_intermediate_files_flag) {
        &purge_files(@to_delete);
    }    
    @to_delete = ();        
    
    
    
    exit(0);
}


####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";

    my $ret = system($cmd);

    if ($ret) {
        confess "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


    
####
sub build_full_paths {
    my ($path, $start_dir) = @_;
    
    my @paths;
    
    foreach my $p (split(/,/, $path)) {

        if ($p !~ /^\//) {
            $p = "$start_dir/$p";
        }
        push (@paths, $p);
    }
    
    $path = join(",", @paths);
    
    return($path);
}

####
sub purge_files {
    my @to_delete = @_;
    
    foreach my $file (@to_delete) {
        if (-e $file || -l $file) {
            print STDERR "-cleaning up and removing intermediate file: $file\n";
            unlink($file);
        }
        else {
            #print STDERR "-warning: cannot locate file $file targeted for deletion.\n";
        }
    }

    return;
}
