#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case bundling);


$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script

my $usage = <<_EOUSAGE_;

################################################################################################################
#
#  --left and --right <string>   (if paired reads)
#     or
#  --single <string>             (if unpaired reads)
#
#  Required inputs:
#   
#  --target <string>            multi-fasta file containing the target sequences (should be named {refName}.fa )
#
#  --seqType <string>           fa | fq    (fastA or fastQ format)
#
# Optional:
#
#  --SS_lib_type <string>       strand-specific library type:  single: F or R  paired: FR or RF
#                                   examples:  single RNA-Ligation method:  F
#                                              single dUTP method: R
#                                              paired dUTP method: RF
#
#  --output|-o <string>                  output directory (default blat_out_dir)
#
# 
#  --num_top_hits <int>         (default: 20) 
#
#  --retain_SAM_files              retain the final SAM files in addition to the BAM files.
#  --retain_intermediate_files     retain all the intermediate sam files produced (they take up lots of space! and there's lots of them)
#  --stop_at_coordSorted.sam.+     stop the run when reaching the coordSorted.sam.+.bam file
#
#  --max_dist_between_pairs             default (2000) 
#
#
#  -I <int>                            maximum intron length  (default: 10000)
#  -P <int>                            min percent identity based on full sequence length  (default: 95)
#  --trim_short_terminal_segments     (trim off short terminal alignment segments that are mostly noise. Default: 10) 
#
#  ## sort-related options
#
#  --sort_buffer_size <string>  (default: 2G)
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
my $max_intron = 10000;
my $min_intron = 20;
my $output_directory = "";
my $trim_short_terminal_segment_length = 10;
my $min_per_ID = 95;
my $num_top_hits = 20;
my $max_dist_between_pairs = 2000;
my $seqType;
my $SS_lib_type;

my $sort_buffer_size = '2G';
my $retain_intermediate_files_flag = 0;
my $retain_SAM_files_flag = 0;
my $NO_BAM = 0;
my $PREP_RSEM = 0;

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
              
              ## Optional:
              "SS_lib_type=s" => \$SS_lib_type,
              
              "I=i" => \$max_intron,
              "i=i" => \$min_intron,
              
              'output|o=s' => \$output_directory,
              
              'trim_short_terminal_segments=i' => \$trim_short_terminal_segment_length,

              'P=i' => \$min_per_ID,
              
              'num_top_hits=i' => \$num_top_hits,
                            
              'max_dist_between_pairs=i' => \$max_dist_between_pairs,
              
              'retain_intermediate_files' => \$retain_intermediate_files_flag,
              'retain_SAM_files' => \$retain_SAM_files_flag,
              
              ## sort options
              'sort_buffer_size=s' => \$sort_buffer_size,
              

              # misc
              'no_bam' => \$NO_BAM,
              'prep_rsem' => \$PREP_RSEM,
    );






if ($help_flag) { die $usage; }

unless ($target_db && -s $target_db) { 
    die $usage . "Must specify target_db and it must exist at that location";
}

unless ($seqType && $seqType =~ /^(fq|fa)$/) {
    die $usage . ", sorry do not understand seqType $seqType";
}


unless ( ($single_file && -e $single_file)
         || 
         ($left_file && -e $left_file
          && $right_file && -e $right_file)) {
    
    if ($single_file) {
        die "Error, cannot find $single_file";
    }
    else {

        die $usage . "sorry, cannot find left: $left_file and right: $right_file";
    }
    
}


unless ($output_directory) {
    $output_directory = "blat_out_dir";
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
    die "Error, SS_lib_type must be one of the following: (F, R, FR, RF)  ";
}

my $UTIL_DIR = "$FindBin::RealBin/../support_scripts";
my $BLAT_UTIL_DIR = "$FindBin::RealBin/blat_util";

my ($start_dir, $work_dir);

main: {
    $start_dir = cwd();

    $single_file = &build_full_paths($single_file, $start_dir) if $single_file;
    $left_file = &build_full_paths($left_file, $start_dir) if $left_file;
    $right_file = &build_full_paths($right_file, $start_dir) if $right_file;
    $target_db = &build_full_paths($target_db, $start_dir);
    
    
    if ($output_directory =~ /^\//) {
        $work_dir = $output_directory;
    }
    else {
        $work_dir = "$start_dir/$output_directory";
    }


    
    &process_cmd("mkdir -p $work_dir") unless (-d $work_dir);
    
    
    my @entries;
    
    if ($single_file) {
        push (@entries, ["single_dir", "single_fa", $single_file]);
    }
    else {
        push (@entries, 
              ["left_dir", "left_fa", $left_file],
              ["right_dir", "right_fa", $right_file]);
    }
    
    chdir $work_dir or die "Error, cannot cd to $work_dir";
    
    unless (-s "target.fa") {
        
        # prep the target here for converting sam to bam later.
    
        my $cmd = "ln -s $target_db target.fa";
        &process_cmd($cmd);
    
        unless (-s "$target_db.fai") {
            my $cmd = "samtools faidx $target_db";
            &process_cmd($cmd);
        }
        $cmd = "ln -s $target_db.fai target.fa.fai";
        &process_cmd($cmd);
    }
    
    ## Run each of the fragment read alignments separately (and simultaneously), and join them later
    
    my @process_monitor_files;
    
    
    foreach my $fa_info (@entries) {
        
        ## always resume work in the work_dir
        chdir $work_dir or die "Error, cannot cd to $work_dir";
        
        my ($target_dir, $trans_fa_name, $trans_file) = @$fa_info;
        
        my $process_monitor_file = &build_full_paths("$target_dir/process_info.txt", $work_dir);
        push (@process_monitor_files, $process_monitor_file);
        
        
        my $pid = fork();
        if ($pid) {
            # parent process
            next;
        }
        
        # child does the work:
        eval {
            unless (-d $target_dir) {
                mkdir $target_dir or die "Error, cannot mkdir $target_dir";
            }
            
            
            # prep target dir with symlinks to components needed.
            if ($seqType eq "fq") {
                my $cmd = "$UTIL_DIR/fastQ_to_fastA.pl -I $trans_file > $target_dir/$trans_fa_name";
                &process_cmd($cmd) unless (-e "$target_dir/$trans_fa_name");
                
            }
            else {
                ## already fasta format or running in bowtie mode, which can use fq format
                my $cmd = "ln -s $trans_file $target_dir/$trans_fa_name";
                &process_cmd($cmd) unless (-e "$target_dir/$trans_fa_name");
            }
            
            
            ## Work in Target_dir
            chdir ($target_dir) or die "Error, cannot cd to $target_dir";
            
            my $cmd = "ln -s ../target.* .";
            &process_cmd($cmd) unless (-e "target.fa");
            
            
            ## processes below ultimately end up generating $target_fa.nameSorted.sam
            
            
            ## run BLAT alignment pipeline
            
            &run_BLAT_alignment_pipeline($trans_fa_name, $trans_file);
            
            
            
            # write monitor result file
            open (my $ofh, ">$process_monitor_file") or die "Error, cannot write to $process_monitor_file";
            
            if ($@) {
                print $ofh "ERROR\t$@\n";
            }
            else {
                print $ofh "SUCCESS\n";
            }
            close $ofh;
            
            exit(0); # child exits. Parent never reaches this statement.
        };
        
    }

    # wait for children to stop running alignments.
    while (wait() > 0) {
        print "-child alignment process completed.\n";
    }
    
    
    chdir $work_dir or die "Error, cannot cd to $work_dir";
    
    ## ensure that each alignment process completed successfully.
    
    my $failure = 0;
    foreach my $monitor_file (@process_monitor_files) {
        if (! -s $monitor_file) {
            print STDERR "Error, cannot find process monitor file: $monitor_file";
            $failure = 1;
            next;
        }
        open (my $fh, $monitor_file);
        my $status_info = <$fh>;
        chomp $status_info;
        my ($status, $msg) = split(/\t/, $status_info);
        if ($status eq "SUCCESS") {
            # no op
        }
        else {
            print STDERR $status_info;
            $failure = 1;
        }
    }
    
    if ($failure) {
        die "Error, alignment pipeline failed due to errors encountered above.";
    }
    else {
        print "\n## Alignment steps succeeded.\n\n";
    }
    
    my @to_delete;
    
    ## merge into single sam file, setting flags properly
    
    my $outfile_basename = basename($output_directory);
    
    if ($single_file) {
        
        my $cmd = "sort -T . -S $sort_buffer_size -k 3,3 -k 4,4n single_dir/single_fa.nameSorted.sam > single_dir/single.coordSorted.sam";
        &process_cmd($cmd) unless (-e "single_dir/single.coordSorted.sam.finished");
        $cmd = "touch single_dir/single.coordSorted.sam.finished";
        &process_cmd($cmd) unless (-e "single_dir/single.coordSorted.sam.finished");
        
        push (@to_delete, "single_dir/single_fa.nameSorted.sam");
        
        $cmd = "cp single_dir/single.coordSorted.sam $outfile_basename.pre.coordSorted.sam";
        &process_cmd($cmd) unless (-e "$outfile_basename.pre.coordSorted.sam.finished");
        $cmd = "touch $outfile_basename.pre.coordSorted.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.pre.coordSorted.sam.finished");
        
        push (@to_delete, "single_dir/single.coordSorted.sam");
        push (@to_delete, "$outfile_basename.pre.coordSorted.sam");
    }
    
    else {
        ## paired mode:
        
        ## Now, capture just the top number of hits taking into account read pairing info.
        
        my $cmd = "$UTIL_DIR/merge_left_right_nameSorted_SAMs.pl --left_sam left_dir/left_fa.nameSorted.sam --right_sam right_dir/right_fa.nameSorted.sam -D $max_dist_between_pairs > combined.nameSorted.sam";
        &process_cmd($cmd) unless (-e "combined.nameSorted.sam.finished");
        $cmd = "touch combined.nameSorted.sam.finished";
        &process_cmd($cmd) unless (-e "combined.nameSorted.sam.finished");
        push (@to_delete, "left_dir/left_fa.nameSorted.sam", "right_dir/right_fa.nameSorted.sam");
        push (@to_delete, "combined.nameSorted.sam");
        
        ## sort by coordinate.
        $cmd = "sort -T . -S $sort_buffer_size -k 3,3 -k 4,4n combined.nameSorted.sam > $outfile_basename.pre.coordSorted.sam";
        &process_cmd($cmd) unless (-e "$outfile_basename.pre.coordSorted.sam.finished");
        $cmd = "touch $outfile_basename.pre.coordSorted.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.pre.coordSorted.sam.finished");
        push (@to_delete, "$outfile_basename.pre.coordSorted.sam");
        
        
    }
    
    
    
            
    
    # report splice junctions and remove short terminal exons that are more likely noise.
    my $cmd = "$FindBin::RealBin/../../Inchworm/bin/cigar_tweaker $outfile_basename.pre.coordSorted.sam target.fa $trim_short_terminal_segment_length | sort -T . -S $sort_buffer_size -k 3,3 -k 4,4n >  $outfile_basename.coordSorted.spliceAdjust.sam";
    &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.spliceAdjust.sam.finished");
    $cmd = "touch $outfile_basename.coordSorted.spliceAdjust.sam.finished";
    &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.spliceAdjust.sam.finished");
    
    push (@to_delete, "$outfile_basename.coordSorted.spliceAdjust.sam");
    
    
    
    # add transcribed orientation info:
    if ($SS_lib_type) {
        my $cmd = "$UTIL_DIR/SAM_set_transcribed_orient_info.pl $outfile_basename.coordSorted.spliceAdjust.sam $SS_lib_type > $outfile_basename.coordSorted.sam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.finished");
        $cmd = "touch $outfile_basename.coordSorted.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.finished");
        
    }
    else {
        # not strand-specific, keep as is and don't disrupt current flow (so use expected output name)
        my $cmd = "cp  $work_dir/$outfile_basename.coordSorted.spliceAdjust.sam $outfile_basename.coordSorted.sam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.finished");
        
        $cmd = "touch $outfile_basename.coordSorted.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.finished");
        
        
    }
    
    ## do some intermediate clean-up
    unless ($retain_intermediate_files_flag) {
        &purge_files(@to_delete);
    }
    @to_delete = (); # reinit        
    
    unless ($NO_BAM) {
        # convert to bam format
        my $cmd = "samtools view -bt target.fa.fai -S $outfile_basename.coordSorted.sam | samtools sort -  $outfile_basename.coordSorted";
        &process_cmd($cmd) if ( (! -e "$outfile_basename.coordSorted.bam.finished") && -s "$outfile_basename.coordSorted.sam");
        $cmd = "touch $outfile_basename.coordSorted.bam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.finished");
        
        push (@to_delete, "$outfile_basename.coordSorted.sam") unless ($retain_SAM_files_flag || $NO_BAM);
    }
    
    ## provide name-sorted SAM
    $cmd = "sort -T . -S $sort_buffer_size -k 1,1 -k 3,3 $outfile_basename.coordSorted.sam > $outfile_basename.nameSorted.sam";
    &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.sam.finished");
    $cmd = "touch $outfile_basename.nameSorted.sam.finished";
    &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.sam.finished");
    
    push (@to_delete, "$outfile_basename.nameSorted.sam") unless ($retain_SAM_files_flag || $NO_BAM);
    
    
    $cmd = "samtools view -bt target.fa.fai $outfile_basename.nameSorted.sam > $outfile_basename.nameSorted.bam";
    &process_cmd($cmd) if  ( (! -e "$outfile_basename.nameSorted.bam") && -s "$outfile_basename.nameSorted.sam" && ! $NO_BAM);
    
    
    $cmd = "samtools index $outfile_basename.coordSorted.bam";
    &process_cmd($cmd) if (-s "$outfile_basename.coordSorted.bam" && (! -s "$outfile_basename.coordSorted.bam.bai") && ! $NO_BAM);
    
    
    if ($SS_lib_type) {
        
        ## strand-specific
        ## separate the sam based on strand, and create separate bam files.  (for convenience sake)
        
        $cmd = "$UTIL_DIR/SAM_strand_separator.pl $outfile_basename.coordSorted.sam $SS_lib_type";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.+.sam.finished" && -e "$outfile_basename.coordSorted.sam.-.sam.finished");
        $cmd = "touch $outfile_basename.coordSorted.sam.+.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.+.sam.finished");
        $cmd = "touch $outfile_basename.coordSorted.sam.-.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.sam.-.sam.finished");
        
        
        $cmd = "$UTIL_DIR/SAM_strand_separator.pl $outfile_basename.nameSorted.sam $SS_lib_type";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.sam.+.sam.finished" && -e "$outfile_basename.nameSorted.sam.-.sam.finished");
        
        $cmd = "touch $outfile_basename.nameSorted.sam.+.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.sam.+.sam.finished");
        $cmd = "touch $outfile_basename.nameSorted.sam.-.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.sam.-.sam.finished");
        
        
        foreach my $sam_file ("$outfile_basename.coordSorted.sam.+.sam", "$outfile_basename.coordSorted.sam.-.sam",
                              "$outfile_basename.nameSorted.sam.+.sam", "$outfile_basename.nameSorted.sam.-.sam") {
            
            push (@to_delete, $sam_file) unless ($retain_SAM_files_flag);
            
            if (-s $sam_file) {
                
                my $bam_file = $sam_file;
                $bam_file =~ s/\.sam$//; # add suffix below
                
                $cmd = ($bam_file =~ /coordSorted/) 
                    ? "samtools view -bt target.fa.fai $sam_file | samtools sort - $bam_file" # .bam ext added auto
                    : "samtools view -bt target.fa.fai $sam_file > $bam_file.bam"; # explicitly adding .bam extension
                
                $bam_file .= ".bam";
                
                &process_cmd($cmd) unless (-e "$bam_file.finished" || $NO_BAM);
                $cmd = "touch $bam_file.finished";
                &process_cmd($cmd) unless (-e "$bam_file.finished" || $NO_BAM);
                
                
                $cmd = "samtools index $bam_file";
                &process_cmd($cmd) if ($bam_file =~ /coordSorted/ && (! -s "$bam_file.bai") && ! $NO_BAM );
                
                
            }
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
sub make_RSEM_bam {
    my ($pre_bam_file) = @_;
    
    
    my $out_bam_prefix = $pre_bam_file;
    $out_bam_prefix =~ s/\.pre\.bam//;

    if ($PREP_RSEM) {
                
        my $cmd = "convert-sam-for-rsem $pre_bam_file $out_bam_prefix -T . ";
        &process_cmd($cmd);
        
        unlink($pre_bam_file);
    }
    else {
        rename($pre_bam_file, $out_bam_prefix) or die "Error, cannot rename $pre_bam_file to $out_bam_prefix";
    }
        
    return;
}

    
    
####
sub build_full_paths {
    my ($path, $start_dir) = @_;
    
    if ($path && $path !~ /^\//) {
        $path = "$start_dir/$path";
    }

    return($path);
}


####
sub run_BLAT_alignment_pipeline {
    my ($trans_fa_name, $trans_file) = @_;
        
    ## Prep sequences>
    if ($seqType eq "fq") {
        my $cmd = "$UTIL_DIR/fastQ_to_tab.pl -I $trans_file > $trans_fa_name.tab";
        &process_cmd($cmd) unless (-s "$trans_fa_name.tab");
    }
    else {
        my $cmd = "$UTIL_DIR/fasta_to_tab.pl  $trans_file NO_FULL_HEADER > $trans_fa_name.tab";
        &process_cmd($cmd) unless (-s "$trans_fa_name.tab");
    }
    
    my $cmd = "sort -T . -S $sort_buffer_size -k 1,1 -k 3,3 $trans_fa_name.tab > $trans_fa_name.sorted.tab";
    &process_cmd($cmd) unless (-s "$trans_fa_name.sorted.tab");
    
    ## run blat
    $cmd = "$BLAT_UTIL_DIR/run_BLAT_shortReads.pl target.fa $trans_fa_name $max_intron $trans_fa_name.psl";
    &process_cmd($cmd) unless (-s "$trans_fa_name.psl");
    
    ## convert to sam
    $cmd = "psl2sam.pl -q 0 -r 0 $trans_fa_name.psl > $trans_fa_name.psl.sam";
    &process_cmd($cmd) unless (-s "$trans_fa_name.psl.sam");
    
    ## sort by name
    $cmd = "sort -T . -S $sort_buffer_size -k 1,1 -k 3,3 $trans_fa_name.psl.sam > $trans_fa_name.psl.nameSorted.sam";
    &process_cmd($cmd) unless (-s "$trans_fa_name.psl.nameSorted.sam");
    
    ## add sequences to sam file.
    $cmd = "$BLAT_UTIL_DIR/blat_sam_add_reads2.pl $trans_fa_name.psl.nameSorted.sam $trans_fa_name.sorted.tab > $trans_fa_name.psl.nameSorted.wReads.sam";
    &process_cmd($cmd) unless (-s "$trans_fa_name.psl.nameSorted.wReads.sam");
    
    ## capture top hits
    $cmd = "$BLAT_UTIL_DIR/top_blat_sam_extractor.pl $trans_fa_name.psl.nameSorted.wReads.sam $num_top_hits $min_per_ID > $trans_fa_name.nameSorted.sam"; # capture twice as many hits in single mode, distill based on pairs later.
    &process_cmd($cmd) unless (-s "$trans_fa_name.nameSorted.sam");
 
    
    return;
    
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
