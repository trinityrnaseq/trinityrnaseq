#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("FindBin::Bin/../PerlLib", "$FindBin::Bin/PerlLib");
use Fasta_reader;
use Getopt::Std;
use strict;
use Carp;
use Cwd;
use Bsub;
use List::Util qw (shuffle);

our ($opt_d, $opt_q, $opt_s, $opt_Q, $opt_b, $opt_p, $opt_O, $opt_h, $opt_c, $opt_B, $opt_X, $opt_M, $opt_S, $opt_o);

&getopts ('dq:s:bp:O:hbc:B:Q:XM:S:o:');

my $usage =  <<_EOH_;

############################# Options ###############################
#
# Required:
# -q query multiFastaFile (full or relative path)
# -p program command line template:   eg. "/path/to/prog [opts] __QUERY_FILE__ [other opts]"
# -o outdir
#
# Optional:
# -c cmds per node 
# -S number of fasta seqs per job submission (default: 1)
# -B bin size  (input seqs per directory)  (default 5000)
# -Q bsub queue
# -M memory   (4000 = 4G is the default setting).  use nubmers only... don't say 4G!!!
# -X commands only!  don't launch.
#
###################### Process Args and Options #####################

_EOH_

    
    ;


if ($opt_h) {
    die $usage;
}

my $CMDS_ONLY = $opt_X;

my $bin_size = $opt_B || 5000;

my $queue = $opt_Q || "week";

our $DEBUG = $opt_d;

my $CMDS_PER_NODE = $opt_c;

my $num_seqs_per_job = $opt_S || 1;

unless ($opt_q && $opt_p && $opt_o) {
    die $usage;
}

my $queryFile = $opt_q;
unless ($queryFile =~ /^\//) {
    $queryFile = cwd() . "/$queryFile";
}

my $program_cmd_template = $opt_p;
unless ($program_cmd_template =~ /__QUERY_FILE__/) {
    die "Error, program cmd template must include '__QUERY_FILE__' placeholder in the command";
}


my $memory = $opt_M || 4000;

my $out_dir = $opt_o;

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

print STDERR "Sequences to search: @searchFileList\n";
my $numFiles = @searchFileList;
print STDERR "There are $numFiles blast search jobs to run.\n";

my $curr_dir = cwd;

if  ($numFiles) {
    
    my @cmds;
    ## formulate blast commands:
    foreach my $searchFile (@searchFileList) {
        $searchFile = "$curr_dir/$searchFile";
        
        my $cmd = $program_cmd_template;
        $cmd =~ s/__QUERY_FILE__/$searchFile/g;
        

        $cmd .= " > $searchFile.OUT ";
		
        unless ($CMDS_ONLY) {
			$cmd .= " 2>$searchFile.ERR";
		}
        push (@cmds, $cmd);
    }
    
	
    @cmds = shuffle(@cmds);

    my $cmds_per_node;
    if ($CMDS_PER_NODE) {
        $cmds_per_node = $CMDS_PER_NODE;
    }
    else {
        $cmds_per_node = int ( scalar(@cmds) / 400); # split job across complete set of nodes available.
        if ($cmds_per_node < 1) {
            # use 10 as default.
            $cmds_per_node = 1;
        }
    }
    

	open (my $fh, ">cmds.list") or die $!;
	foreach my $cmd (@cmds) {
		print $fh "$cmd\n";
	}
	close $fh;
	
	unless ($CMDS_ONLY) {
		my $bsubber = new Bsub({cmds=>\@cmds,
								log_dir => $out_dir,
								cmds_per_node => $cmds_per_node,
								queue => $queue,
								memory => $memory,
							}
			);
		
		$bsubber->bsub_jobs();
		
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
