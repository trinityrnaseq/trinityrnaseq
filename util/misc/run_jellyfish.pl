#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use FindBin;

my $usage = "\n\nusage: $0 reads.fa hash_size\n\n";

my $reads_file = $ARGV[0] or die $usage;
my $hash_size = $ARGV[1] or die $usage;


my $JELLYFISH_DIR = $FindBin::RealBin . "/../../trinity-plugins/jellyfish-1.1.3";
my $CPU = 4;
my $min_kmer_cov = 1;

unless ($reads_file =~ /^\//) {
    $reads_file = cwd() . "/$reads_file";
}


my $workdir = "H_" . ($hash_size/1e9) . "G";
mkdir($workdir) or die "Error, cannot mkdir $workdir";
chdir ($workdir) or die "Error, cannot cd to $workdir";


my $jelly_kmer_fa_file = "jellyfish.kmers.fa";
            
#    my $jelly_hash_size = int( ($max_memory - $read_file_size)/7); # decided upon by Rick Westerman
             
my $cmd = "$JELLYFISH_DIR/bin/jellyfish count -t $CPU -m 25 -s $hash_size ";

#        $cmd .= " --both-strands ";
            
$cmd .= " $reads_file";
            
&process_cmd($cmd);

my @kmer_db_files;

foreach my $file (<mer_counts_*>) {
    my $cmd = "$JELLYFISH_DIR/bin/jellyfish dump -L $min_kmer_cov $file >> $file.kmer_fa";
    &process_cmd($cmd);
    
    $cmd = "cat $file.kmer_fa >> $jelly_kmer_fa_file";
    &process_cmd($cmd);
    
}
            
$cmd = "grep '>' $jelly_kmer_fa_file | wc -l | tee kmer_count.txt";
&process_cmd($cmd);


exit(0);


####
sub process_cmd {
    my ($cmd) = @_;

    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


