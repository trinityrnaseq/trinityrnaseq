#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;

my $usage = "usage: $0 fasta\n";

my $fasta_file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);


my %node_id_to_seq;

while (my $seq_obj = $fasta_reader->next()) {
	
	my $header = $seq_obj->get_header();
	my $seq = $seq_obj->get_sequence();
    my $accession = $seq_obj->get_accession();

    my $gene_id = $accession;
    $gene_id =~ s/_i\d+$//;
    

    my $path_info;

    if ($header =~ /path=\[([^\]]+)\]/) {
        $path_info = $1;
    }
    else {
        die "Error, cannot extract path info from $header";
    }

	#print "path: $path_info\n";
    
    my @nodes = split(/\s+/, $path_info);
    for my $node (@nodes) {
        my ($node_id, $seq_range) = split(/:/, $node);
        
        $node_id = join(":", $gene_id, $node_id);
        
        my ($lend, $rend) = split(/-/, $seq_range);
        my $node_seq = substr($seq, $lend, $rend - $lend + 1);

        #print "$node_id\t$node_seq\n";
        
        if (exists $node_id_to_seq{$node_id}) {
            if ($node_id_to_seq{$node_id} ne $node_seq) {
                
                my ($longer_seq, $shorter_seq) = reverse sort {length($a)<=>length($b)} ($node_seq, $node_id_to_seq{$node_id});
                if (index($longer_seq, $shorter_seq) > 0) {
                    # keep the longer one.
                    $node_id_to_seq{$node_id} = $longer_seq;
                }
                else {
                    print STDERR "-warning, difference in node seqs: $node_id\t$node_id_to_seq{$node_id}\tvs\t$node_seq\n";
                }
            }
        }
        else {
            $node_id_to_seq{$node_id} = $node_seq;
        }
        
    }
}


for my $node_id (sort keys %node_id_to_seq) {
    my $node_seq = $node_id_to_seq{$node_id};
    print "$node_id\t" . length($node_seq) . "\t$node_seq\n";
}


exit(0);

