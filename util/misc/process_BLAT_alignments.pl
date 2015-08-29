#!/usr/bin/env perl

use FindBin;
use strict;
use warnings;
use threads;
use POSIX qw(ceil);
use lib "$FindBin::Bin/../../PerlLib";

use Fasta_reader;
use Process_cmd;
use Thread_helper;

use Getopt::Long qw (:config no_ignore_case bundling);
use vars qw ($DEBUG $opt_h $opt_g $opt_t $opt_c $opt_o $opt_B $opt_I);

my $CPU = 1;
my $num_top_hits = 1;
my $KEEP_PSLX = 0;

&GetOptions( 'g=s' => \$opt_g,
             'd' => \$DEBUG,
             'h' => \$opt_h,
             'c=s' => \$opt_c,
             'o=s' => \$opt_o,
             't=s' => \$opt_t,
			 'I=i' => \$opt_I,
             'N=i' => \$num_top_hits,
             'CPU=i' => \$CPU,
             'KEEP_PSLX' => \$KEEP_PSLX,
             );

our $SEE = 0;

$|++;

my $MAX_INTRON = $opt_I || 100000;

my $usage = <<_EOH_;

Script chunks EST alignments into more manageable data sets.

############################# Options ###############################
#
# -g <string>       genomic_seq.db
# -t <string>       transcripts database
# -I <int>          maximum intron length (default: 100000)
# -o <string>       prefix for output file (default: 'blat')
# -N <int>          number of top hits (default: $num_top_hits)
#
# --CPU <int>       number of threads (default: 1)
#
# -h this help menu
# -d debug mode
# --KEEP_PSLX      retain the raw blat output files
#
###################### Process Args and Options #####################


_EOH_

    ;


my $genome_db = $opt_g;
my $transcript_db = $opt_t;
my $output_prefix = $opt_o || "blat";
my $blat_path = "blat";
my $util_dir = $FindBin::Bin;

unless ($genome_db && $transcript_db) {
    die "$usage\n";
}

my $ooc_cmd = "$blat_path $genome_db $transcript_db -q=rna -dots=100 -maxIntron=$MAX_INTRON  -makeOoc=11.ooc /dev/null";

my $blat_thr;

unless (-s "11.ooc") {
    $blat_thr = threads->create('process_cmd', $ooc_cmd);
}

my $blat_partitions_dir = "blat_out_dir";
unless (-d $blat_partitions_dir) {
    mkdir($blat_partitions_dir) or die "Error, cannot mkdir $blat_partitions_dir";
}

my $num_seqs = `grep '>' $transcript_db | wc -l `;
$num_seqs =~ s/\s//g;
unless ($num_seqs && $num_seqs =~ /^\d+$/) {
    die "Error, cannot determine number of fasta entries in $transcript_db";
}

my $seqs_per_partition = ceil($num_seqs/$CPU);
unless ($seqs_per_partition > 1) { 
    $seqs_per_partition = 1;
}

my @transcript_files = &partition_transcript_db($transcript_db, $seqs_per_partition, $blat_partitions_dir);

if ($blat_thr) {
    $blat_thr->join();
    if ($blat_thr->error()) {
        die "Error, $ooc_cmd died ...";
    }
}


###############################
## Run BLAT
###############################
my @pslx_files;

{

    my $thread_helper = new Thread_helper($CPU);
    my %thread_to_checkpoint;
    
    foreach my $transcript_file (@transcript_files) {
        
        
        ## process blat search:
        my $cmd = "$blat_path $genome_db $transcript_file -q=rna -dots=100 "
            . " -maxIntron=$MAX_INTRON -out=pslx -ooc=11.ooc $transcript_file.pslx";
        
        
        my $checkpoint_file = "$transcript_file.pslx.completed";
        unless (-e $checkpoint_file) {
            
            my $thread = threads->create('process_cmd', $cmd);
            $thread_helper->add_thread($thread);
            
            my $thread_id = $thread->tid();
            $thread_to_checkpoint{$thread_id} = {thread => $thread,
                                                 checkpoint => $checkpoint_file,
                                             };
        }
        
        push (@pslx_files, "$transcript_file.pslx");
        
    }    
    
    $thread_helper->wait_for_all_threads_to_complete();

    ## write checkpoints for successful threads:
    foreach my $thread_id (keys %thread_to_checkpoint) {
        my $thread_info_href = $thread_to_checkpoint{$thread_id};
        my $thread = $thread_info_href->{thread};
        unless ($thread->error()) {
            my $checkpoint = $thread_info_href->{checkpoint};
            system("touch $checkpoint");
        }
    }
    
    if (my @failures = $thread_helper->get_failed_threads()) {
        die "Error, ". scalar(@failures) . " blat searches failed. ";
    }
}

######################
## get top hits.
######################


my @top_hits_files;


{
    
    my %thread_to_checkpoint;
    my $thread_helper = new Thread_helper($CPU);

    foreach my $pslx_file (@pslx_files) {
        
        my $cmd = "$util_dir/blat_top_hit_extractor.pl $pslx_file $num_top_hits > $pslx_file.top_${num_top_hits}";
        
        my $completed_checkpoint_file = "$pslx_file.top_${num_top_hits}.completed";
        unless (-e $completed_checkpoint_file) {
            
            my $thread = threads->create('process_cmd', $cmd);
            $thread_helper->add_thread($thread);
            
            my $thread_id = $thread->tid();
            $thread_to_checkpoint{$thread_id} = { thread => $thread,
                                                  checkpoint => $completed_checkpoint_file,
                                              };
            
        }
        push (@top_hits_files, "$pslx_file.top_${num_top_hits}");
    }
    
    $thread_helper->wait_for_all_threads_to_complete();
    
    ## write checkpoints for successful threads:
    foreach my $thread_id (keys %thread_to_checkpoint) {
        my $thread_info_href = $thread_to_checkpoint{$thread_id};
        my $thread = $thread_info_href->{thread};
        unless ($thread->error()) {
            my $checkpoint = $thread_info_href->{checkpoint};
            system("touch $checkpoint");
        }
    }
    

    if (my @failures = $thread_helper->get_failed_threads()) {
        die "Error, ". scalar(@failures) . " blat top hit selectors failed. ";
    }
}



if (-s "$output_prefix.gff3") {
    print STDERR "WARNING, REPLACING EXISTING FILE: $output_prefix.gff3.  KILL THIS NOW TO PREVENT THIS. (you have 10 seconds)\n";
    sleep(10);
    print STDERR "OK, too late.  replacing it now.\n";
    unlink("$output_prefix.gff3");
}

foreach my $top_hits_file (@top_hits_files) {
    
    # convert to gff3 format
    print STDERR "-converting $top_hits_file to gff3\n";
        
    my $cmd = "$util_dir/pslx_to_gff3.pl < $top_hits_file >> $output_prefix.gff3";
    &process_cmd($cmd);
}

## clean up the pslx files we no longer need.
foreach my $pslx_file (@pslx_files) {
    unlink($pslx_file) unless $KEEP_PSLX; # these files can be huge. Once have top hits, no longer need all hits (hopefully).
}

print STDERR "done.\n";

exit(0);


####
sub partition_transcript_db {
    my ($transcript_db, $seqs_per_partition, $blat_partitions_dir) = @_;
    
    my @files;
    
    my $checkpoint_file = "$blat_partitions_dir/partitions.completed";

    if (-e $checkpoint_file) {
        my @files = <$blat_partitions_dir/partition.*.fa>;
        return(@files);
    }
    else {

        
        my $fasta_reader = new Fasta_reader($transcript_db);
        my $partition_counter = 0;
        
        my $counter = 0;
        my $ofh;
        
        while (my $seq_obj = $fasta_reader->next()) {
            my $fasta_entry = $seq_obj->get_FASTA_format();
            if ($counter % $seqs_per_partition == 0) {
                close $ofh if $ofh;
                $partition_counter++;
                my $outfile = "$blat_partitions_dir/partition.$counter.fa";
                open ($ofh, ">$outfile") or die "Error, cannot write to outfile: $outfile";
                push (@files, $outfile);
            }
            print $ofh $fasta_entry;
            $counter++;
        }
        
        close $ofh if $ofh;
        
        system("touch $checkpoint_file");
        
        return(@files);
    }
}
