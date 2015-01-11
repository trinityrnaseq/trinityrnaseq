#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "\n\n\tusage: $0 fasta_file min_seq_len (byGene|byTrans) [numPerDir=100]\n\nnote: byGene restricts to multi-iso genes\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $min_seq_length = $ARGV[1] or die $usage;
my $by_gene_or_trans = $ARGV[2] or die $usage;
my $num_per_dir = $ARGV[3] || 100;

unless ($by_gene_or_trans =~ /^(byGene|byTrans)$/) { die $usage; }

main: {
    
    my $fasta_reader = new Fasta_reader($fasta_file);
    my %seqs = $fasta_reader->retrieve_all_seqs_hash();
    
    my $by_gene_flag = ($by_gene_or_trans eq 'byGene');

    %seqs = &reorganize_seqs(\%seqs, $by_gene_flag);

    my $outdir = "$by_gene_or_trans.dir";

    mkdir($outdir) or die "Error, cannot mkdir $outdir";
    
    open (my $ofh_file_listing, ">$outdir.listing") or die "Error, cannot write to $outdir.listing";
    
    my $counter = 0;
    foreach my $acc (keys %seqs) {
        my $bindir = "$outdir/bin." . int($counter/$num_per_dir);
        if (! -d $bindir) {
            mkdir ($bindir) or die "Error, cannot mkdir $bindir";
        }
        
        my $acc_for_filename = $acc;
        $acc_for_filename =~ s/\W/_/g;

        my @seq_entries = @{$seqs{$acc}};

        if ($by_gene_flag && scalar(@seq_entries) == 1) { 
            ## only focusing on multi-iso genes.
            next;
        }


        my $filename = "$bindir/$acc_for_filename.ref.fa";
        open (my $ofh, ">$filename") or die "Error, cannot write to $filename";
        

        foreach my $seq_entry (@seq_entries) {
            my ($acc, $seq) = @$seq_entry;
            print $ofh ">$acc\n$seq\n";
        }
        close $ofh;
        
        print $ofh_file_listing "$filename\n";
        print STDERR "// wrote $filename\n";
        
        $counter++;
    }
    
    
    print STDERR "\n\nDone.\n\n";
    
    exit(0);
    

}



####
sub reorganize_seqs {
    my ($fasta_seqs_href, $by_gene_flag) = @_;
    
    my %reorg_fasta;
    
    foreach my $acc (sort keys %$fasta_seqs_href) {
        
        my $seq = uc $fasta_seqs_href->{$acc};
        
        $seq =~ s/N//g;  # no N characters allowed.
        
        unless (length($seq) >= $min_seq_length) { next; }
        
        my $key = $acc;
        if ($by_gene_flag) {
            my ($trans, $gene) = split(/;/, $acc);
            unless ($gene) {
                confess "Error, no gene ID extracted from $acc ";
            }
            $key = $gene;
        }
        
                
        push (@{$reorg_fasta{$key}}, [$acc, $seq]);
    }

    return(%reorg_fasta);
    
}

    
