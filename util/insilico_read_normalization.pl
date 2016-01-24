#!/usr/bin/env perl

use strict;
use warnings;
use threads;
no strict qw(subs refs);

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use File::Basename;
use Cwd;
use Carp;
use Getopt::Long qw(:config no_ignore_case pass_through);
use Fastq_reader;
use Fasta_reader;
use threads;
use Data::Dumper;
use COMMON;

open (STDERR, ">&STDOUT");  ## capturing stderr and stdout in a single stdout stream


## Jellyfish
my $max_memory;

# Note: For the Trinity logo below the backslashes are quoted in order to keep
#   them from quoting the character than follows them.  "\\" keeps "\ " from occuring.

my $output_directory = cwd();
my $help_flag;
my $seqType;
my @left_files;
my @right_files;
my $left_list_file;
my $right_list_file;
my @single_files;
my $SS_lib_type;
my $CPU = 2;
my $MIN_KMER_COV_CONST = 2;  ## DO NOT CHANGE
my $max_cov;
my $pairs_together_flag = 0;
my $max_pct_stdev = 200;
my $KMER_SIZE = 25;

my $__devel_report_kmer_cov_stats = 0;

my $PARALLEL_STATS = 0;
my $JELLY_S;

my $usage = <<_EOUSAGE_;


###############################################################################
#
# Required:
#
#  --seqType <string>      :type of reads: ( 'fq' or 'fa')
#  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for 
#                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
#                     
#
#  --max_cov <int>         :targeted maximum coverage for reads.
#
#
#  If paired reads:
#      --left  <string>    :left reads   (if specifying multiple files, list them as comma-delimited. eg. leftA.fq,leftB.fq,...)
#      --right <string>    :right reads
#
#  Or, if unpaired reads:
#      --single <string>   :single reads
#
#  Or, if you have read collections in different files you can use 'list' files, where each line in a list
#  file is the full path to an input file.  This saves you the time of combining them just so you can pass
#  a single file for each direction.
#      --left_list  <string> :left reads, one file path per line
#      --right_list <string> :right reads, one file path per line
#
####################################
##  Misc:  #########################
#
#  --pairs_together                :process paired reads by averaging stats between pairs and retaining linking info.
#
#  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
#                                   if paired: RF or FR,
#                                   if single: F or R.   (dUTP method = RF)
#                                   See web documentation.
#  --output <string>               :name of directory for output (will be
#                                   created if it doesn't already exist)
#                                   default( "${output_directory}" )
#
#  --CPU <int>                     :number of threads to use (default: = $CPU)
#  --PARALLEL_STATS                :generate read stats in parallel for paired reads
#
#  --KMER_SIZE <int>               :default $KMER_SIZE
#
#  --max_pct_stdev <int>           :maximum pct of mean for stdev of kmer coverage across read (default: $max_pct_stdev)
#
#  --no_cleanup                    :leave intermediate files                      
#  --tmp_dir_name <string>         default("tmp_normalized_reads");
#
###############################################################################




_EOUSAGE_

    ;

my $ROOTDIR = "$FindBin::RealBin/../";
my $UTILDIR = "$ROOTDIR/util/support_scripts/";
my $INCHWORM_DIR = "$ROOTDIR/Inchworm";
my $JELLYFISH_DIR = "$ROOTDIR/trinity-plugins/jellyfish";
my $FASTOOL_DIR = "$ROOTDIR/trinity-plugins/fastool";

unless (@ARGV) {
    die "$usage\n";
}

my $NO_FASTOOL = 0;

my $NO_CLEANUP = 0;

my $TMP_DIR_NAME = "tmp_normalized_reads";


&GetOptions( 
             
    'h|help' => \$help_flag,

    ## general opts
    "seqType=s" => \$seqType,
    "left=s{,}" => \@left_files,
    "right=s{,}" => \@right_files,
    "single=s{,}" => \@single_files,
    
    "left_list=s" => \$left_list_file,
    "right_list=s" => \$right_list_file,
    
    "SS_lib_type=s" => \$SS_lib_type,
    "max_cov=i" => \$max_cov,
    "output=s" => \$output_directory,

    # Jellyfish
    'JM=s'          => \$max_memory, # in GB

    # misc
    'no_fastool' => \$NO_FASTOOL,

    'KMER_SIZE=i' => \$KMER_SIZE,
    'CPU=i' => \$CPU,
    'PARALLEL_STATS' => \$PARALLEL_STATS,
    'kmer_size=i' => \$KMER_SIZE,
    'max_pct_stdev=i' => \$max_pct_stdev,
    'pairs_together' => \$pairs_together_flag,

     'no_cleanup' => \$NO_CLEANUP,

     #devel
     '__devel_report_kmer_cov_stats' => \$__devel_report_kmer_cov_stats,
    'jelly_s=i' => \$JELLY_S,

    'tmp_dir_name=s' => \$TMP_DIR_NAME, 
    
);



if ($help_flag) {
    die "$usage\n";
}

if (@ARGV) {
    die "Error, do not understand options: @ARGV\n";
}


my $USE_FASTOOL = 1; # by default, using fastool for fastq to fasta conversion
if ($NO_FASTOOL) {
    $USE_FASTOOL = 0;
}


unless ($seqType =~ /^(fq|fa)$/) {
    die "Error, set --seqType to 'fq' or 'fa'";
}
unless ($max_memory && $max_memory =~ /^\d+G/) {
    die "Error, must set --JM to number of G of RAM (ie. 10G)  ";
}

if ($SS_lib_type) {
    unless ($SS_lib_type =~ /^(R|F|RF|FR)$/) {
        die "Error, unrecognized SS_lib_type value of $SS_lib_type. Should be: F, R, RF, or FR\n";
    }

    ## note, if single-end reads and strand-specific, just treat it as F and don't bother revcomplementing here (waste of time)
    if ($SS_lib_type eq 'R') {
        $SS_lib_type = 'F';
    }
    
}


if ($left_list_file) {
    @left_files = &read_list_file($left_list_file);
}
if ($right_list_file) {
    @right_files = &read_list_file($right_list_file);
}

unless (@single_files || (@left_files && @right_files)) {
    die "Error, need either options 'left' and 'right' or option 'single'\n";
}


if (@left_files) {
    @left_files = split(",", join(",", @left_files));
}
if (@right_files) {
    @right_files = split(",", join(",", @right_files));
}
if (@single_files) {
    @single_files = split(",", join(",", @single_files));
}



unless ($max_cov && $max_cov >= 2) {
    die "Error, need to set --max_cov at least 2";
}



## keep the original 'xG' format string for the --JM option, then calculate the numerical value for max_memory
my $JM_string = $max_memory;    ## this one is used in the Chrysalis exec string
my $sort_mem;
if ($max_memory) {
    $max_memory =~ /^([\d\.]+)G$/ or die "Error, cannot parse max_memory value of $max_memory.  Set it to 'xG' where x is a numerical value\n";
    


    $max_memory = $1;
    
    # prep the sort memory usage 
    $sort_mem = $max_memory;
    if ($PARALLEL_STATS) {
        $sort_mem = int($sort_mem/2);
    }
    $sort_mem .= "G";
    
    $max_memory *= 1024**3; # convert to from gig to bytes
}
else {
    die "Error, must specify max memory for jellyfish to use, eg.  --JM 10G \n";
}

if ($pairs_together_flag && ! ( @left_files && @right_files) ) {
    die "Error, if setting --pairs_together, must use the --left and --right parameters.";
}


my $sort_exec = &COMMON::get_sort_exec($CPU);



main: {
    
    my $start_dir = cwd();

    ## create complete paths for input files:
    @left_files = &create_full_path(@left_files) if @left_files;
    @right_files = &create_full_path(@right_files) if @right_files;
    $left_list_file = &create_full_path($left_list_file) if $left_list_file;
    $right_list_file = &create_full_path($right_list_file) if $right_list_file;
    @single_files = &create_full_path(@single_files) if @single_files;
    $output_directory = &create_full_path($output_directory);
    
    unless (-d $output_directory) {
        
        mkdir $output_directory or die "Error, cannot mkdir $output_directory";
    }
    
    chdir ($output_directory) or die "Error, cannot cd to $output_directory";

    my $tmp_directory = "$output_directory/$TMP_DIR_NAME";
    
    my $CREATED_TMP_DIR_HERE_FLAG = 0;
    if (! -d $tmp_directory) {
        mkdir $tmp_directory or die "Error, cannot mkdir $tmp_directory";
        $CREATED_TMP_DIR_HERE_FLAG = 1;
    }
    chdir $tmp_directory or die "Error, cannot cd to $tmp_directory";
    
    my $trinity_target_fa = (@single_files) ? "single.fa" : "both.fa"; 
    
    my @files_need_stats;
    my @checkpoints;


    if ( (@left_files && @right_files) || 
         ($left_list_file && $right_list_file) ) {
        
        my ($left_SS_type, $right_SS_type);
        if ($SS_lib_type) {
            ($left_SS_type, $right_SS_type) = split(//, $SS_lib_type);
        }

        print("Converting input files. (both directions in parallel)");

        my $thr1;
        my $thr2;

                
        if (-s "left.fa" && -e "left.fa.ok") {
            $thr1 = threads->create(sub { print ("left file exists, nothing to do");});
        }
        else {
            $thr1 = threads->create('prep_list_of_seqs', \@left_files, $seqType, "left", $left_SS_type);
            push (@checkpoints, ["left.fa", "left.fa.ok"]);
        }
        
        if (-s "right.fa" && -e "right.fa.ok") {
            $thr2 = threads->create(sub { print ("right file exists, nothing to do");});
        }
        else {
            $thr2 = threads->create('prep_list_of_seqs', \@right_files, $seqType, "right", $right_SS_type);
            push (@checkpoints, ["right.fa", "right.fa.ok"]);
        }
		
        $thr1->join();
		$thr2->join();

        if ($thr1->error() || $thr2->error()) {
            die "Error, conversion thread failed";
        }
        
        &process_checkpoints(@checkpoints);
        
		print("Done converting input files.");
        
        push (@files_need_stats, 
              [\@left_files, "left.fa"], 
              [\@right_files, "right.fa"]);
        
        @checkpoints = ();
        &process_cmd("cat left.fa right.fa > $trinity_target_fa") unless (-s $trinity_target_fa && -e "$trinity_target_fa.ok");
        unless (-s $trinity_target_fa == ((-s "left.fa") + (-s "right.fa"))){
            die "$trinity_target_fa (".(-s $trinity_target_fa)." bytes) is different from the combined size of left.fa and right.fa (".((-s "left.fa") + (-s "right.fa"))." bytes)\n";
        }
        push (@checkpoints, [ $trinity_target_fa, "$trinity_target_fa.ok" ]);
        
        
    } 
    elsif (@single_files) {

        ## Single-mode

        unless (-s "single.fa" && -e "single.fa.ok") {
            &prep_list_of_seqs(\@single_files, $seqType, "single", $SS_lib_type);

        }
        push (@files_need_stats, [\@single_files, "single.fa"]);
        push (@checkpoints, [ "single.fa", "single.fa.ok" ]);
        
    } 
    else {
        die "not sure what to do. "; # should never get here.
    }

    &process_checkpoints(@checkpoints);
    
    my $kmer_file = &run_jellyfish($trinity_target_fa, $SS_lib_type);

    &generate_stats_files(\@files_need_stats, $kmer_file, $SS_lib_type);

    if ($pairs_together_flag) {
        &run_nkbc_pairs_together(\@files_need_stats, $kmer_file, $SS_lib_type, $max_cov, $max_pct_stdev);
    } else {
        &run_nkbc_pairs_separate(\@files_need_stats, $kmer_file, $SS_lib_type, $max_cov, $max_pct_stdev);
    }
        
    my @outputs;
    @checkpoints = ();
    my @threads;
    foreach my $info_aref (@files_need_stats) {
        my ($orig_file, $converted_file, $stats_file, $selected_entries) = @$info_aref;

        ## do multi-threading

        my $base = (scalar @$orig_file == 1) ? basename($orig_file->[0]) : basename($orig_file->[0]) . "_ext_all_reads";
        
        my $normalized_filename_prefix = $output_directory . "/$base.normalized_K${KMER_SIZE}_C${max_cov}_pctSD${max_pct_stdev}";
        my $outfile;
        
        if ($seqType eq 'fq') {
            $outfile = "$normalized_filename_prefix.fq";
        }
        else {
            # fastA
            $outfile = "$normalized_filename_prefix.fa";
        }
        push (@outputs, $outfile);
    
        ## run in parallel
        
        my $checkpoint_file = "$outfile.ok";
        unless (-e $checkpoint_file) {
            my $thread = threads->create('make_normalized_reads_file', $orig_file, $seqType, $selected_entries, $outfile);
                
            push (@threads, $thread);
            push (@checkpoints, [$outfile, $checkpoint_file]);
        }
    }
    
    my $num_fail = 0;
    foreach my $thread (@threads) {
        $thread->join();
        if ($thread->error()) {
            print STDERR "Error encountered with thread.\n";
            $num_fail++;
        }
    }
    if ($num_fail) {
        die "Error, at least one thread died";
    }
    
    &process_checkpoints(@checkpoints);
    
    chdir $output_directory or die "Error, cannot chdir to $output_directory";
    
    ## link them up with simpler names so they're easy to find by downstream scripts.
    if (scalar @outputs == 2) {
        my $left_out = $outputs[0];
        my $right_out = $outputs[1];
        
        &process_cmd("ln -sf $left_out left.norm.$seqType");
        &process_cmd("ln -sf $right_out right.norm.$seqType");
    }
    else {
        my $single_out = $outputs[0];
        &process_cmd("ln -sf $single_out single.norm.$seqType");
    }
    
    unless ($NO_CLEANUP) {
        if ($CREATED_TMP_DIR_HERE_FLAG) {
            print STDERR "-removing tmp dir $tmp_directory\n";
            `rm -rf $tmp_directory`;
        }
    }
    
    print "\n\nNormalization complete. See outputs: \n\t" . join("\n\t", @outputs) . "\n";
    

    exit(0);
}


####
sub build_selected_index {
    my $file = shift;
    
    my %index = ();
    
    open(my $ifh, $file) || die "failed to read selected_entries file $file: $!";
    
    while (my $line = <$ifh> ) {
        chomp $line;
        next unless $line =~ /\S/;
        
        ## want core, .... just in case.
        $line =~ s|/\w$||;
        
        #print STDERR "-want $line\n";


        $index{$line} = 0;
    }
    
    return \%index;
}


####
sub make_normalized_reads_file {
    my ($source_files_aref, $seq_type, $selected_entries, $outfile) = @_;

    open (my $ofh, ">$outfile") or die "Error, cannot write to $outfile";

    my @source_files = @$source_files_aref;
    
    my $idx = build_selected_index( $selected_entries );
    
    for my $orig_file ( @source_files ) {
        my $reader;
        
        # if we had a consistent interface for the readers, we wouldn't have to code this up separately below... oh well.
        ##  ^^ I enjoyed this lamentation, so I left it in the rewrite - JO
        if    ($seqType eq 'fq') { $reader = new Fastq_reader($orig_file) } 
        elsif ($seqType eq 'fa') { $reader = new Fasta_reader($orig_file) }
        else {  die "Error, do not recognize format: $seqType" }
        
        while ( my $seq_obj = $reader->next() ) {
        
            my $acc;
        
            if ($seqType eq 'fq') {
                $acc = $seq_obj->get_core_read_name();
            } elsif ($seqType eq 'fa') {
                $acc = $seq_obj->get_accession();
                $acc =~ s|/[12]$||;
            }
            
            if ( exists $$idx{$acc} ) {
                $$idx{$acc}++;
                my $record = '';
                
                if    ($seqType eq 'fq') { $record = $seq_obj->get_fastq_record() } 
                elsif ($seqType eq 'fa') { $record = $seq_obj->get_FASTA_format(fasta_line_len => -1) }
                
                print $ofh $record;
            }
        }
    }
    
    ## check and make sure they were all found
    my $not_found_count = 0;
    for my $k ( keys %$idx ) {
        $not_found_count++ if $$idx{$k} == 0;
    }
    
    if ( $not_found_count ) {
        die "Error, not all specified records have been retrieved (missing $not_found_count) from @source_files";
    }
    
    return;
}


####
sub run_jellyfish {
    my ($reads, $strand_specific_flag) = @_;
    
    my $jelly_kmer_fa_file = "jellyfish.K${KMER_SIZE}.min${MIN_KMER_COV_CONST}.kmers.fa";
    
    print STDERR "-------------------------------------------\n"
        . "----------- Jellyfish  --------------------\n"
        . "-- (building a k-mer catalog from reads) --\n"
        . "-------------------------------------------\n\n";
    
    my $jellyfish_checkpoint = "$jelly_kmer_fa_file.success";
    
    unless (-e $jellyfish_checkpoint) {


        my $read_file_size = -s $reads;
        
        my $jelly_hash_size = int( ($max_memory - $read_file_size)/7); # decided upon by Rick Westerman
        
        
        if ($jelly_hash_size < 100e6) {
            $jelly_hash_size = 100e6; # seems reasonable for a min hash size as 100M
        }

        ## for testing
        if ($JELLY_S) {
            $jelly_hash_size = $JELLY_S;
        }
                
        my $cmd = "$JELLYFISH_DIR/bin/jellyfish count -t $CPU -m $KMER_SIZE -s $jelly_hash_size ";
        
        unless ($SS_lib_type) {
            ## count both strands
            $cmd .= " --canonical ";
        }
        
        $cmd .= " $reads";
        
        &process_cmd($cmd);
        
           
        if (-s $jelly_kmer_fa_file) {
            unlink($jelly_kmer_fa_file) or die "Error, cannot unlink $jelly_kmer_fa_file";
        }

        my $jelly_db = "mer_counts.jf";
        
        ## write a histogram of the kmer counts.
        $cmd = "$JELLYFISH_DIR/bin/jellyfish histo -t $CPU -o $jelly_kmer_fa_file.histo $jelly_db";
        &process_cmd($cmd);



        $cmd = "$JELLYFISH_DIR/bin/jellyfish dump -L $MIN_KMER_COV_CONST $jelly_db > $jelly_kmer_fa_file";

        &process_cmd($cmd);
        
        unlink($jelly_db);
            
        ## if got this far, consider jellyfish done.
        &process_cmd("touch $jellyfish_checkpoint");
        
    }
    

    return($jelly_kmer_fa_file);
}


####  (from Trinity.pl)
## WARNING: this function appends to the target output file, so a -s check is advised
#   before you call this for the first time within any given script.
sub prep_seqs {
    my ($initial_file, $seqType, $file_prefix, $SS_lib_type) = @_;

    ($initial_file) = &add_fifo_for_gzip($initial_file) if $initial_file =~ /\.gz$/;
    
    if ($seqType eq "fq") {
        # make fasta
        
        my $perlcmd = "$UTILDIR/fastQ_to_fastA.pl -I $initial_file ";
        my $fastool_cmd = "$FASTOOL_DIR/fastool";
        if ($SS_lib_type && $SS_lib_type eq "R") {
            $perlcmd .= " --rev ";
            $fastool_cmd .= " --rev ";
        }
        $fastool_cmd .= " --illumina-trinity --to-fasta $initial_file >> $file_prefix.fa";
        $perlcmd .= " >> $file_prefix.fa";  
        
       
        my $cmd = ($USE_FASTOOL) ? $fastool_cmd : $perlcmd;
        
        &process_cmd($cmd);
    }
    elsif ($seqType eq "fa") {
        if ($SS_lib_type && $SS_lib_type eq "R") {
            my $cmd = "$UTILDIR/revcomp_fasta.pl $initial_file >> $file_prefix.fa";
            &process_cmd($cmd);
        }
        else {
            ## just symlink it here:
            my $cmd = "ln -s $initial_file $file_prefix.fa";
            &process_cmd($cmd) unless (-e "$file_prefix.fa");
        }
    }
    elsif (($seqType eq "cfa") | ($seqType eq "cfq")) {
        # make double-encoded fasta
        my $cmd = "$UTILDIR/csfastX_to_defastA.pl -I $initial_file ";
        if ($SS_lib_type && $SS_lib_type eq "R") {
            $cmd .= " --rev ";
        }
        $cmd .= ">> $file_prefix.fa";
        &process_cmd($cmd);
  }
    return;
}



###
sub prep_list_of_seqs {
    my ($files, $seqType, $file_prefix, $SS_lib_type) = @_;


    for my $file ( @$files ) {
        prep_seqs( $file,  $seqType, $file_prefix, $SS_lib_type);
    }
    
    return 0;
}


###
sub create_full_path {
    my (@files) = @_;

    my @ret;
    foreach my $file (@files) {
        my $cwd = cwd();
        if ($file !~ m|^/|) { # must be a relative path
            $file = $cwd . "/$file";
        }
        push (@ret, $file);
    }
    
    if (wantarray) {
        return(@ret);
    }
    else {
        if (scalar @ret > 1) {
            confess("Error, provided multiple files as input, but only requesting one file in return");
        }
        return($ret[0]);
    }
}


###
sub read_list_file {
    my ($file, $regex) = @_;
    
    my @files;
    
    open(my $ifh, $file) || die "failed to read input list file ($file): $!";
    
    while (my $line = <$ifh>) {
        chomp $line;
        next unless $line =~ /\S/;
        
        if ( defined $regex ) {
            if ( $line =~ /$regex/ ) {
                push @files, $line;
            }
        } else {
            push @files, $line;
        }
    }
    
    return @files;
}


####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $start_time = time();
    my $ret = system("bash", "-c", $cmd);
    my $end_time = time();

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    print "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

    return;
}

####
sub generate_stats_files {
    my ($files_need_stats_aref, $kmer_file, $SS_lib_type) = @_;
    
    my @cmds;

    my @checkpoints;
    foreach my $info_aref (@$files_need_stats_aref) {
        my ($orig_file, $converted_fa_file) = @$info_aref;

        my $stats_filename = "$converted_fa_file.K$KMER_SIZE.stats";
        push (@$info_aref, $stats_filename);
        
        my $cmd = "$INCHWORM_DIR/bin/fastaToKmerCoverageStats --reads $converted_fa_file --kmers $kmer_file --kmer_size $KMER_SIZE  --num_threads $CPU ";
        unless ($SS_lib_type) {
            $cmd .= " --DS ";
        }
        
        if ($__devel_report_kmer_cov_stats) {
            $cmd .= " --capture_coverage_info ";
        }
                
        $cmd .= " > $stats_filename";
    
        push (@cmds, $cmd) unless (-e "$stats_filename.ok");
        push (@checkpoints, ["$stats_filename", "$stats_filename.ok"]);
    }
    
    if (@cmds) {
        if ($PARALLEL_STATS) {
            &process_cmds_parallel(@cmds);
        }
        else {
            &process_cmds_serial(@cmds);
        }
    }
    
    &process_checkpoints(@checkpoints);
    

    
    {
        ## sort by read name
        print STDERR "-sorting each stats file by read name.\n";
        my @cmds;
        @checkpoints = ();
        foreach my $info_aref (@$files_need_stats_aref) {
            my $stats_file = $info_aref->[-1];
            my $sorted_stats_file = $stats_file . ".sort";
            my $cmd = "$sort_exec -k5,5 -T . -S $sort_mem $stats_file > $sorted_stats_file";
            push (@cmds, $cmd) unless (-e "$sorted_stats_file.ok");
            $info_aref->[-1] = $sorted_stats_file;
            push (@checkpoints, [$sorted_stats_file, "$sorted_stats_file.ok"]);
            
        }
        
        if (@cmds) {
            if ($PARALLEL_STATS) {
                &process_cmds_parallel(@cmds);
            }
            else {
                &process_cmds_serial(@cmds);
            }
            
        }

        &process_checkpoints(@checkpoints);        
        
    }
    
    return;
}


####
sub process_checkpoints {
    my @checkpoints = @_;
    
    foreach my $checkpoint (@checkpoints) {
        my ($outfile, $checkpoint_file) = @$checkpoint;
        if (-s "$outfile" && ! -e $checkpoint_file) {
            &process_cmd("touch $checkpoint_file");
        }
    }
    
    return;
}


####
sub run_nkbc_pairs_separate {
    my ($files_need_stats_aref, $kmer_file, $SS_lib_type, $max_cov, $max_pct_stdev) = @_;

    my @cmds;

    my @checkpoints;
    foreach my $info_aref (@$files_need_stats_aref) {
        my ($orig_file, $converted_file, $stats_file) = @$info_aref;
                
        my $selected_entries = "$stats_file.C$max_cov.pctSD$max_pct_stdev.accs";
        my $cmd = "$UTILDIR/nbkc_normalize.pl $stats_file $max_cov $max_pct_stdev > $selected_entries";
        push (@cmds, $cmd) unless (-e "$selected_entries.ok");

        push (@$info_aref, $selected_entries);
        
        push (@checkpoints, [$selected_entries, "$selected_entries.ok"]);

    }


    &process_cmds_parallel(@cmds); ## low memory, all I/O - fine to always run in parallel.

    &process_checkpoints(@checkpoints);
    
    return;
        
}


####
sub run_nkbc_pairs_together {
    my ($files_need_stats_aref, $kmer_file, $SS_lib_type, $max_cov, $max_pct_stdev) = @_;

    my $left_stats_file = $files_need_stats_aref->[0]->[2];
    my $right_stats_file = $files_need_stats_aref->[1]->[2];
        
    my $pair_out_stats_filename = "pairs.K$KMER_SIZE.stats";
    
    my $cmd = "$UTILDIR/nbkc_merge_left_right_stats.pl --left $left_stats_file --right $right_stats_file --sorted";

    $cmd .= " > $pair_out_stats_filename";
    
    &process_cmd($cmd) unless (-e "$pair_out_stats_filename.ok");
    my @checkpoints = ( [$pair_out_stats_filename, "$pair_out_stats_filename.ok"] );
    &process_checkpoints(@checkpoints);
    

    unless (-s $pair_out_stats_filename) {
        die "Error, $pair_out_stats_filename is empty.  Be sure to check your fastq reads and ensure that the read names are identical except for the /1 or /2 designation.";
    }
    
    my $selected_entries = "$pair_out_stats_filename.C$max_cov.pctSD$max_pct_stdev.accs";
    $cmd = "$UTILDIR/nbkc_normalize.pl $pair_out_stats_filename $max_cov $max_pct_stdev > $selected_entries";
    &process_cmd($cmd) unless (-e "$selected_entries.ok");
    @checkpoints = ( [$selected_entries, "$selected_entries.ok"] );
    &process_checkpoints(@checkpoints);

    push (@{$files_need_stats_aref->[0]}, $selected_entries);
    push (@{$files_need_stats_aref->[1]}, $selected_entries);
    
    
    return;
    
}



####
sub process_cmds_parallel {
    my @cmds = @_;


    my @threads;
    foreach my $cmd (@cmds) {
        # should only be 2 cmds max
        my $thread = threads->create('process_cmd', $cmd);
        push (@threads, $thread);
    }
                
    my $ret = 0;
    
    foreach my $thread (@threads) {
        $thread->join();
        if (my $error = $thread->error()) {
            print STDERR "Error, thread exited with error $error\n";
            $ret++;
        }
    }
    if ($ret) {
        die "Error, $ret threads errored out";
    }

    return;
}

####
sub process_cmds_serial {
    my @cmds = @_;

    foreach my $cmd (@cmds) {
        &process_cmd($cmd);
    }

    return;
}


    
####
sub add_fifo_for_gzip {
    my @files = @_;

    foreach my $file (@files) {
        if (ref $file eq "ARRAY") {
            my @f = &add_fifo_for_gzip(@{$file});
            $file = [@f];
        }
        elsif ($file =~ /\.gz$/) {
            $file = "<(gunzip -c $file)";
        }
    }
    
    return(@files);

}
