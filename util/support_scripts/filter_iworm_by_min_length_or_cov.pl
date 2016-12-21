#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;


my $usage = "usage: $0 IwormFastaFile min_length min_cov\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $min_length = $ARGV[1] or die $usage;
my $min_cov = $ARGV[2] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

while (my $seqobj = $fasta_reader->next()) {
    my $fasta_entry = $seqobj->get_FASTA_format();
    my $sequence = $seqobj->get_sequence();

    my $iworm_acc = $seqobj->get_accession();
    my ($iworm_num, $iworm_cov, @rest) = split(/;/, $iworm_acc);
    
    unless (length($sequence) >= $min_length || $iworm_cov >= $min_cov) {
        next;
    }
    
    print $fasta_entry;
}


exit(0);

