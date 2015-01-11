#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use Cwd;
use File::Basename;
use Carp;
use Data::Dumper;

use Getopt::Long qw(:config no_ignore_case bundling);

$ENV{PATH} .= "\:$FindBin::Bin/../trinity-plugins/rsem/sam/";  # include samtools in path, already included in rsem build.

$ENV{LC_ALL} = 'C';  # critical for proper sorting using [system "sort -k1,1 ..."] within the perl script

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
#  --aligner <string>           bowtie, bowtie2
#
# Optional:
#
#  --SS_lib_type <string>       strand-specific library type:  single: F or R  paired: FR or RF
#                                   examples:  single RNA-Ligation method:  F
#                                              single dUTP method: R
#                                              paired dUTP method: RF
#
#  --output|-o <string>                  output directory (default \${aligner}_out)
#
# 
#  --num_top_hits <int>         (default: 20) 
#
#  --retain_intermediate_files     retain all the intermediate sam files produced (they take up lots of space! and there's lots of them)
#  --prep_rsem                     prep the rsem-ready files
#  --run_rsem                      execute rsem (implies --prep_rsem)
#        --trinity_mode       extract gene/trans mapping info from Trinity.fasta file directly
#        --gene_trans_map <string>    rsem gene-to-transcript mapping file to use.
#
#     --max_dist_between_pairs             default (2000) 
#
#  --just_prep_build               just prepare the bowtie-build and stop.
#
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
#     $0 --seqType fq --left left.fq --right right.fq --target Trinity.fasta --aligner bowtie -- -p 4 --all --best --strata -m 300
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

my $output_directory = "";
my $num_top_hits = 20;
my $max_dist_between_pairs = 2000;
my $seqType;
my $SS_lib_type;
my $aligner;
my $retain_intermediate_files_flag = 0;
my $PREP_RSEM = 0;
my $RUN_RSEM = 0;

my $trinity_mode;
my $gene_trans_map_file;

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
              
              "aligner=s" => \$aligner,


              ## Optional:
              "SS_lib_type=s" => \$SS_lib_type,
              
              'output|o=s' => \$output_directory,
              
              'num_top_hits=i' => \$num_top_hits,
                            
              'max_dist_between_pairs=i' => \$max_dist_between_pairs,
              
              'retain_intermediate_files' => \$retain_intermediate_files_flag,
              
              'prep_rsem' => \$PREP_RSEM,
              'run_rsem' => \$RUN_RSEM,
              'trinity_mode' => \$trinity_mode,
              'gene_trans_map=s' =>\$gene_trans_map_file,
              
              
              'just_prep_build' => \$JUST_PREP_BUILD,
              
              );


if ($help_flag) { die $usage; }

if ($RUN_RSEM) {
    $PREP_RSEM = 1;
}


unless ($target_db && -s $target_db) { 
    die $usage . "Must specify target_db and it must exist at that location";
}

unless ($JUST_PREP_BUILD || ($seqType && $seqType =~ /^(fq|fa)$/)) {
    die $usage . ", sorry do not understand seqType $seqType";
}

unless ($aligner && $aligner =~ /^bowtie|bowtie2$/) {
    die $usage . "sorry, aligner [$aligner] is not supported.";
}

unless ($JUST_PREP_BUILD 
        || 
        ($left_file && $right_file)
        ||
        ($single_file)

    )  {
    die $usage . "sorry, cannot find files $left_file $right_file $single_file";
}


unless ($output_directory) {
    $output_directory = "${aligner}_out";
}

if ($SS_lib_type && $SS_lib_type !~ /^(F|R|FR|RF)$/) {
    die "Error, SS_lib_type must be one of the following: (F, R, FR, RF)  ";
}


unless ($aligner eq "bowtie") {
    
    if ($PREP_RSEM) {
        die "Error, can only prep RSEM files if aligner = bowtie";
    }
}



## check for required programs
{
    
    my @required_progs = qw(samtools);
    
    if ($aligner eq "bowtie" || $aligner eq "tophat") {
        push (@required_progs, qw (bowtie-build bowtie) );
    }
    elsif ($aligner eq "bowtie2" || $aligner eq "tophat2") {
        push (@required_progs, qw (bowtie2-build bowtie2));
    }
    

    if ($PREP_RSEM) {
        push (@required_progs, "convert-sam-for-rsem");
    }
    

    foreach my $prog (@required_progs) {
        my $path = `which $prog`;
        unless ($path =~ /^\//) {
            die "Error, path to required $prog cannot be found";
        }
    }
}


my $util_dir = "$FindBin::Bin/../util/support_scripts";

my ($start_dir, $work_dir, $num_hits);


main: {
    $start_dir = cwd();

    $left_file = &build_full_paths($left_file, $start_dir) if $left_file;
    $right_file = &build_full_paths($right_file, $start_dir) if $right_file;
    $target_db = &build_full_paths($target_db, $start_dir);
    
    $single_file = &build_full_paths($single_file, $start_dir) if $single_file;
    
    
    if ($output_directory =~ /^\//) {
        $work_dir = $output_directory;
    }
    else {
        $work_dir = "$start_dir/$output_directory";
    }

    
    &process_cmd("mkdir -p $work_dir") unless (-d $work_dir);


    ######################################
    ## Prep the bowtie index of the target
    ######################################
    
    my $index_ext = ($aligner eq 'bowtie') ? "ebwt" : "bt2";
        
    my @bowtie_build_files = <$target_db.*.$index_ext>;
    print Dumper(\@bowtie_build_files);
    unless (@bowtie_build_files) {
        
        print STDERR "Note - bowtie-build indices do not yet exist. Indexing genome now.\n";
        ## run bowtie-build:
        my $curr_dir = cwd();
        my $target_dir = dirname($target_db);
        my $basename = basename($target_db);
        chdir($target_dir) or die "Error, cannot cd to $target_dir";
        
        my $start_index_file = "$target_db._${index_ext}_idx_.start";
        if (-e $start_index_file) {
            die "Error, indexing in progress or never finished. Wait for other process to finish building the index, or delete the $start_index_file and try again";
        }
        &process_cmd("touch $start_index_file");
        # make bowtie index
        
        my $builder = ($aligner =~ /2$/) ? "bowtie2-build" : "bowtie-build";
        
        my $cmd = "$builder -q $target_db $target_db";
        &process_cmd($cmd);
        
        unlink($start_index_file);
        
        chdir($curr_dir) or die "Error, cannot cd back to $curr_dir";
    }
    
    if ($JUST_PREP_BUILD) {
        print STDERR "Just preparing build, stopping now.\n";
        exit(0);
    }
    
    
    my @entries;
    if ($left_file && $right_file) {
        @entries = (["left", $left_file],
                    ["right", $right_file]);
    }
    else {
        @entries = (["single", $single_file]);
    }
    
    chdir $work_dir or die "Error, cannot cd to $work_dir";
    
    unless (-s "target.fa") {
        
        # prep the target here for converting sam to bam later.
        
    
        my $cmd = "ln -sf $target_db target.fa";
        &process_cmd($cmd);
    
    
        unless (-s "$target_db.fai") {
            my $cmd = "samtools faidx $target_db";
            &process_cmd($cmd);
        }
        $cmd = "ln -sf $target_db.fai target.fa.fai";
        &process_cmd($cmd);
        
            
        @bowtie_build_files = <$target_db.*.$index_ext>;
        
        print Dumper(\@bowtie_build_files);
        
        ## reuse them, but name them target.fa
        foreach my $file (@bowtie_build_files) {
            if ($file =~ /offrate|TRANS/) { next; } # ignore the offrate-1 index

            my @parts = split(/\./, $file);
            my $ebwt_ext = pop @parts;
            my $other_ext = pop @parts;
            my $rev = pop @parts;
            if ($rev eq "rev") {
                $rev = "rev.";
            }
            else {
                $rev = "";
            }
            my $new_filename = "target." . "${rev}" . "${other_ext}.$ebwt_ext";
            my $cmd = "ln -sf $file $new_filename";
            &process_cmd($cmd);
        }
        
        
    }

        
    
    ## Run each of the fragment read alignments separately (and simultaneously), and join them later
    
    $num_hits = 2 * $num_top_hits; 
    
    my @process_monitor_files;
    
    
    foreach my $fa_info (@entries) {
        
        ## always resume work in the work_dir
        chdir $work_dir or die "Error, cannot cd to $work_dir";
        
        my ($target_dir, $trans_file) = @$fa_info;
                

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
            
            
            ## Work in Target_dir
            chdir ($target_dir) or die "Error, cannot cd to $target_dir";
            
            my $cmd = "ln -sf ../target.* .";
            &process_cmd($cmd) unless (-e "target.fa");
            
            
            ## processes below ultimately end up generating $target.nameSorted.sam
            
            &run_bowtie_alignment_pipeline($target_dir, $trans_file);
            
            
        };
                    
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
    
    if ($left_file && $right_file) {
        ## Now, capture just the top number of hits taking into account read pairing info.
        
        my $cmd = "$util_dir/merge_left_right_nameSorted_SAMs.pl --left_sam left/left.nameSorted.bam --right_sam right/right.nameSorted.bam -D $max_dist_between_pairs | samtools view -bt target.fa.fai -S - > combined.nameSorted.pre.bam";
        &process_cmd($cmd) unless (-e "combined.nameSorted.pre.bam.finished");
        $cmd = "touch combined.nameSorted.pre.bam.finished";
        &process_cmd($cmd) unless (-e "combined.nameSorted.pre.bam.finished");
        push (@to_delete, "left_dir/left_fa.nameSorted.bam", "right_dir/right_fa.nameSorted.bam");
        push (@to_delete, "combined.nameSorted.pre.bam");
        
        ## sort by coordinate.
        $cmd = "samtools sort -o combined.nameSorted.pre.bam - > $outfile_basename.coordSorted.pre.bam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.pre.bam.finished");
        $cmd = "touch $outfile_basename.coordSorted.pre.bam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.pre.bam.finished");
        push (@to_delete, "$outfile_basename.coordSorted.pre.bam");
    }
    else {
        ## single file
        my $cmd = "samtools sort -o single/single.nameSorted.bam - > $outfile_basename.coordSorted.pre.bam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.pre.bam.finished");
        $cmd = "touch $outfile_basename.coordSorted.pre.bam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.pre.bam.finished");
        push (@to_delete, "$outfile_basename.coordSorted.pre.bam");
    }
    
    # add transcribed orientation info:
    if ($SS_lib_type) {
        my $cmd = "$util_dir/SAM_set_transcribed_orient_info.pl $outfile_basename.coordSorted.pre.bam $SS_lib_type | samtools view -bt target.fa.fai -S -o - - | samtools sort -o - - >  $outfile_basename.coordSorted.bam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.finished");
        $cmd = "touch $outfile_basename.coordSorted.bam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.finished");
        
    }
    else {
        # not strand-specific, keep as is and don't disrupt current flow (so use expected output name)
        my $cmd = "samtools sort -o $outfile_basename.coordSorted.pre.bam - > $outfile_basename.coordSorted.bam";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.finished");
        
        $cmd = "touch $outfile_basename.coordSorted.bam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.finished");
                
    }
    
    ## do some intermediate clean-up
    unless ($retain_intermediate_files_flag) {
        &purge_files(@to_delete);
    }
    @to_delete = (); # reinit        
    
    ## provide name-sorted SAM
    my $cmd = "samtools sort -no $outfile_basename.coordSorted.bam - > $outfile_basename.nameSorted.bam";
    &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.bam.finished");
    $cmd = "touch $outfile_basename.nameSorted.bam.finished";
    &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.bam.finished");
        
    
    $cmd = "samtools view -bt target.fa.fai $outfile_basename.nameSorted.sam > $outfile_basename.nameSorted.bam";
    &process_cmd($cmd) if  ( (! -e "$outfile_basename.nameSorted.bam") && -s "$outfile_basename.nameSorted.sam");
    
    
    $cmd = "samtools index $outfile_basename.coordSorted.bam";
    &process_cmd($cmd) if (-s "$outfile_basename.coordSorted.bam" && (! -s "$outfile_basename.coordSorted.bam.bai") );


    my @rsem_prep_files;
    
    if ($SS_lib_type) {
        
        ## strand-specific
        ## separate the sam based on strand, and create separate bam files.  (for convenience sake)
        
        $cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.coordSorted.bam $SS_lib_type";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.+.sam.finished" && -e "$outfile_basename.coordSorted.bam.-.sam.finished");
        $cmd = "touch $outfile_basename.coordSorted.bam.+.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.+.sam.finished");
        $cmd = "touch $outfile_basename.coordSorted.bam.-.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.coordSorted.bam.-.sam.finished");
        
        
        $cmd = "$util_dir/SAM_strand_separator.pl $outfile_basename.nameSorted.bam $SS_lib_type";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.bam.+.sam.finished" && -e "$outfile_basename.nameSorted.bam.-.sam.finished");
        
        $cmd = "touch $outfile_basename.nameSorted.bam.+.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.bam.+.sam.finished");
        $cmd = "touch $outfile_basename.nameSorted.bam.-.sam.finished";
        &process_cmd($cmd) unless (-e "$outfile_basename.nameSorted.bam.-.sam.finished");
        
        
        foreach my $sam_file ("$outfile_basename.coordSorted.bam.+.sam", "$outfile_basename.coordSorted.bam.-.sam",
                                  "$outfile_basename.nameSorted.bam.+.sam", "$outfile_basename.nameSorted.bam.-.sam") {
            

            my $bam_file = $sam_file;
            $bam_file =~ s/\.sam$/\.bam/; # add suffix below

            if ($PREP_RSEM && $bam_file =~ /nameSorted\.bam\.\+/) { # only doing the sense mappings
                push (@rsem_prep_files, $bam_file);
            }
            
            unless (-s $sam_file) { next; } # empty file
            
    
            push (@to_delete, $sam_file);
        
            $cmd = ($bam_file =~ /coordSorted/) 
                ? "samtools view -bt target.fa.fai $sam_file | samtools sort -o - - > $bam_file" # .bam ext added auto
                : "samtools view -bt target.fa.fai $sam_file > $bam_file"; # explicitly adding .bam extension
            
            
            &process_cmd($cmd) unless (-e "$bam_file.finished");
            $cmd = "touch $bam_file.finished";
            &process_cmd($cmd) unless (-e "$bam_file.finished");
            
            
            $cmd = "samtools index $bam_file";
            &process_cmd($cmd) if ($bam_file =~ /coordSorted/ && (! -s "$bam_file.bai") );
            
                       
        }
    }
    else { 
        
        ## not strand specific
        
        if ($PREP_RSEM) {
            
            # this only makes sense if paired reads are involved
            
            push (@rsem_prep_files, "$outfile_basename.nameSorted.bam");
        }
        
    }

    
    if ($PREP_RSEM) {
        ## create files for RSEM
        foreach my $bam_file (@rsem_prep_files) {
            
            my $rsem_pre_bam;
            if ($left_file && $right_file) {
                $rsem_pre_bam = "$bam_file.PropMapPairsForRSEM.pre.bam";
                $cmd = "$util_dir/SAM_extract_properly_mapped_pairs.pl $bam_file | samtools view -bt target.fa.fai -S - > $rsem_pre_bam";
                &process_cmd($cmd) unless (-e "$rsem_pre_bam.finished");
                $cmd = "touch $rsem_pre_bam.finished";
                &process_cmd($cmd) unless (-e "$rsem_pre_bam.finished");
                
                push (@to_delete, $rsem_pre_bam, "$rsem_pre_bam.finished");

            }
            else {
                # no prop pair extraction required.
                $rsem_pre_bam = $bam_file;
            }
            
            
            &make_RSEM_bam($rsem_pre_bam) if (-s $rsem_pre_bam);
            # above makes .bam w/o .pre
            

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

    my $rsem_bam = "$out_bam_prefix.bam";
    
    if ($PREP_RSEM) {
                
        my $cmd = "convert-sam-for-rsem $pre_bam_file $out_bam_prefix -T . ";
        &process_cmd($cmd);
        
        
    }
    else {
        $rsem_bam = $pre_bam_file;
        
    }
    
    if ($RUN_RSEM) {
        
        my $cmd = "$FindBin::Bin/align_and_estimate_abundance.pl --est_method RSEM --aln_method $rsem_bam "
            . " --transcripts $target_db --seqType $seqType ";
        
        if ($left_file && $right_file) {
            $cmd .= " --left $left_file --right $right_file ";
        }
        else {
            $cmd .= " --single $single_file ";
        }
        if ($SS_lib_type) {
            $cmd .= " --SS_lib_type $SS_lib_type ";
        }
        if ($trinity_mode) {
            $cmd .= " --trinity_mode ";
        }
        elsif ($gene_trans_map_file) {
            $cmd .= " --gene_trans_map $gene_trans_map_file ";
        }
        &process_cmd($cmd);
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
sub run_bowtie_alignment_pipeline {
    my ($target, $trans_fa_name) = @_;
        
    my $format = ($seqType eq "fq") ? "-q" : "-f";
    
    my @bowtie_opts = @ARGV;
    if ($aligner eq "bowtie" && ! grep { $_ eq "-a" || $_ eq "--all" } @bowtie_opts) {
        push (@bowtie_opts, "-k", $num_hits);  ## don't specifiy the number of hits if using the -a (all) option.
    }
    
    ## FIXME: need to account for use of multiple read file inputs here.
    
    my @adj_files;
    foreach my $file (split(/,/, $trans_fa_name)) {
        
        if ($file =~ /\.gz$/) {
            $file = "<(zcat $file)";
        }
        push (@adj_files, $file);
    }
    
    $trans_fa_name = join(",", @adj_files);
    
    my $cmd = "bash -c \"$aligner @bowtie_opts --chunkmbs 512 -S $format target $trans_fa_name | samtools view -S -b -o $target.pre.bam - \"";

    &process_cmd($cmd) unless (-e "$target.pre.bam.finished");
    $cmd = "touch $target.pre.bam.finished";
    &process_cmd($cmd) unless (-e "$target.pre.bam.finished");
        
    
    ## remove unaligned reads
    $cmd = "samtools view -F 4 -b $target.pre.bam > $target.bam";
    &process_cmd($cmd) unless (-e "$target.bam.finished");
    $cmd = "touch $target.bam.finished";
    &process_cmd($cmd) unless (-e "$target.bam.finished");
    unless ($retain_intermediate_files_flag) {
        &purge_files("$target.pre.bam");
    }
    
    # name-sort 
    $cmd = "samtools sort -no $target.bam - > $target.nameSorted.bam";
    &process_cmd($cmd) unless (-e "$target.nameSorted.bam.finished");
    $cmd = "touch $target.nameSorted.bam.finished";
    &process_cmd($cmd) unless (-e "$target.nameSorted.bam.finished");
    unless ($retain_intermediate_files_flag) {
        &purge_files("$target.bam");
    }
    

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
