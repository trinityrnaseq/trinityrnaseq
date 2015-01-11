#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use strict;
use Carp;
use Cwd;
use HTC::GridRunner;
use List::Util qw (shuffle);
use File::Basename;

my $usage =  <<_EOH_;

############################# Options ###############################
#
##  Required:
#
# --grid_conf|G <string>   grid config file
#
# --query|Q <string>       query multiFastaFile (full or relative path)
# --search|S <string>      search multiFastaFile (full or relative path)
# --program|P <string>     program (e.g. blastn blastx tblastn)
#
##  Optional
#
# --options|O <string>     blast options, default "-max_target_seqs 1 -outfmt 6 -evalue 1e-5"
# 
# --outdir|o <string>      outdir (default: "htc_blast_outdir")
# --fasta_per_job|F <int>  number of fasta seqs per job submission (default: 100)
#
# -X                       commands only!  don't launch.
#
# --RESUME                 resume from earlier search attempt
#
###################### Process Args and Options #####################



_EOH_

    
    ;

my $help_flag;
my $grid_conf_file;
my $queryFile;
my $searchDB;
my $program;
my $progOptions = "-max_target_seqs 1 -outfmt 6 -evalue 1e-5";
my $num_seqs_per_job = 100;
my $out_dir = "htc_blast_outdir";
my $CMDS_ONLY;
our $DEBUG;
my $bin_size = 2000; # files per directory 

my $RESUME_FLAG = 0;

&GetOptions ( 'help|h' => \$help_flag,
              'grid_conf|G=s' => \$grid_conf_file,
              'queryFile|Q=s' => \$queryFile,
              'searchDB|S=s' => \$searchDB,
              'program|P=s' => \$program,
              'options|O=s' => \$progOptions,
              'fasta_per_job|F=i' => \$num_seqs_per_job,
              'out_dir|o=s' => \$out_dir,
              'X' => \$CMDS_ONLY,
              'd' => \$DEBUG,
              'RESUME' => \$RESUME_FLAG,

);



if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, not recognizing parameters: @ARGV ";
}


unless ($queryFile && $searchDB && $program) {
    print STDERR "Need to specify -Q, -S, and -P ";
    die $usage;
}

unless ($queryFile =~ /^\//) {
    $queryFile = cwd() . "/$queryFile";
}

unless ($searchDB =~ /^\//) {
    $searchDB = cwd() . "/$searchDB";
}

unless (-s $searchDB || "-s $searchDB.pal") {
    die "Error, can't find $searchDB\n";
}


unless ($out_dir =~ /^\//) {
    # create full path
    $out_dir = cwd() . "/$out_dir";
}


my $cache_file_prefix = "htc-" . join("-", basename($queryFile), basename($searchDB), $program);
my $cache_file = "$cache_file_prefix.cache_success";
my $cache_cmds = "$cache_file_prefix.cmds";
my @cmds;

if ($RESUME_FLAG) {
    @cmds = `cat $cache_cmds`;
    unless (@cmds) {
        die "Error, cannot resume, no cmds found - expecting file: $cache_cmds";
    }
    chomp @cmds;
}
else {

    ## Create files to search
    
    my $fastaReader = new Fasta_reader($queryFile);
    
    my @searchFileList;
    
    my $count = 0;
    my $current_bin = 1;
    
    mkdir $out_dir or die "Error, cannot mkdir $out_dir";
    
    my $bindir = "$out_dir/grp_" . sprintf ("%04d", $current_bin);
    mkdir ($bindir) or die "Error, cannot mkdir $bindir";
    
    
    
    while (my $fastaSet = &get_next_fasta_entries($fastaReader, $num_seqs_per_job) ) {
        
        $count++;
        
        my $filename = "$bindir/$count.fa";
        
        push (@searchFileList, $filename);
        
        open (TMP, ">$filename") or die "Can't create file ($filename)\n";
        print TMP $fastaSet;
        close TMP;
        chmod (0666, $filename);
    	
        if ($count % $bin_size == 0) {
            # make a new bin:
            $current_bin++;
            $bindir = "$out_dir/grp_" . sprintf ("%04d", $current_bin);
            mkdir ($bindir) or die "Error, cannot mkdir $bindir";
        }
    }
    
    print "Sequences to search: @searchFileList\n";
    my $numFiles = @searchFileList;
    print "There are $numFiles blast search jobs to run.\n";
    
    if  ($numFiles) {
        
        ## formulate blast commands:
        foreach my $searchFile (@searchFileList) {
            
            my $cmd = "$program -db $searchDB -query $searchFile $progOptions > $searchFile.$program.result ";
            unless ($CMDS_ONLY) {
                $cmd .= "2>$searchFile.$program.stderr";
            }
            push (@cmds, $cmd);
        }
        
        open (my $fh, ">$cache_cmds") or die $!;
        foreach my $cmd (@cmds) {
            print $fh "$cmd\n";
        }
        close $fh;
    }
}

if (@cmds) {
    unless ($CMDS_ONLY) {
        
        my $grid_runner = new HTC::GridRunner($grid_conf_file, $cache_file);
        my $ret = $grid_runner->run_on_grid(@cmds);
        
        if ($ret) {
            die "Error, not all butterfly commands could complete successfully... cannot continue.";
        }
    }
    
} else {
    print STDERR "Sorry, no searches to perform.  Results already exist here\n";
}

## Cleanup

exit(0);


####
sub get_next_fasta_entries {
    my ($fastaReader, $num_seqs) = @_;


    my $fasta_entries_txt = "";
    
    for (1..$num_seqs) {
        my $seq_obj = $fastaReader->next();
        unless ($seq_obj) {
            last;
        }

        my $entry_txt = $seq_obj->get_FASTA_format();
        $fasta_entries_txt .= $entry_txt;
    }

    return($fasta_entries_txt);
}
