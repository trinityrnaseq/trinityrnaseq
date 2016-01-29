#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::RealBin/../../PerlLib/");
use Ascii_genome_illustrator;
use Cwd;

my $usage = "usage: $0 refseqs.fa target.fa min_per_id [min_pct_aligned]\n\n";

my $REFSEQS_FA = $ARGV[0] or die $usage;
my $TARGET_FA = $ARGV[1] or die $usage;
my $MIN_PER_ID = $ARGV[2] or die $usage;
my $MIN_PCT_ALIGNED = $ARGV[3] || 0;

main: {

    
    my $cmd = "blat $TARGET_FA $REFSEQS_FA -t=dna -q=dna -out=blast9 $TARGET_FA.blat.out";

    &process_cmd($cmd);

    # tried blastn, but it doesnt like the Trinity.fasta header format...
    #my $cmd = "makeblastdb -in trinity_out_dir/Trinity.fasta -dbtype nucl";    
    #$cmd = " blastn -db trinity_out_dir/Trinity.fasta -query refseqs.fa -outfmt 6 > blast.blastn.out";
    #&process_cmd($cmd);
    
    &generate_ascii_illustration("$TARGET_FA.blat.out");
    
    exit(0);
}

####
sub generate_ascii_illustration {
    my ($align_out_file) = @_;

    my %ref_seq_lengths = &get_seq_lengths($REFSEQS_FA);
    my %target_seq_lengths = &get_seq_lengths($TARGET_FA);
    
    my @hits = &parse_align_out($align_out_file);
    
    ## Generate a reference view.
    foreach my $ref_acc (keys %ref_seq_lengths) {
        my $length = $ref_seq_lengths{$ref_acc};
        my $ascii_illustration = new Ascii_genome_illustrator($ref_acc, 60);
        
        #$ascii_illustration->add_feature($ref_acc, 1, $length, "=");
        
        my @matches = grep { $_->{ref_acc} eq $ref_acc } @hits;
        @matches = reverse sort {$a->{bitscore}<=>$b->{bitscore}} @matches;


        foreach my $match (@matches) {
            
            my $target_acc = $match->{target_acc};
            my $target_start = $match->{target_start};
            my $target_end = $match->{target_end};
            my $per_id = $match->{per_id};
            
            if ($per_id < $MIN_PER_ID) { next; }

            my $ref_start = $match->{ref_start};
            my $ref_end = $match->{ref_end};

            my $target_length = $target_seq_lengths{$target_acc};
            my $pct_of_target_aligned = sprintf("%.2f", (abs($target_end-$target_start) + 1) / $target_length * 100);
            
            if ($pct_of_target_aligned < $MIN_PCT_ALIGNED) { next; }

            my $target_feature_name = "$target_acc $target_start-$target_end:$target_length ($pct_of_target_aligned\% aln, $per_id\% ID)";
            
            $ascii_illustration->add_feature($target_feature_name, $ref_start, $ref_end, "-");
            
        }
        
        print $ascii_illustration->illustrate(1, $length);
        print "\n"; # spacer
        
    }
    
    return;
}

####
sub parse_align_out {
    my ($align_out_file) = @_;
    
    my @alignments;

    open (my $fh, $align_out_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;

=blast_outfmt6

# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score


0       NM_012009_Sh2d1b1
1       comp0_c0_seq1
2       100.00
3       1433
4       0
5       0
6       1
7       1433
8       1
9       1433
10      0.0
11      2647

=cut

        my @x = split(/\t/);
        my ($ref_acc, $target_acc, $per_id, $align_len, $mismatches, $gaps, $ref_start, $ref_end, $target_start, $target_end, $evalue, $bitscore) = @x;
        
        my $struct = { target_acc => $target_acc,
                       ref_acc => $ref_acc,
                       
                       per_id => $per_id,
                       mismatches => $mismatches,
                       gaps => $gaps,
                       
                       target_start => $target_start,
                       target_end => $target_end,

                       ref_start => $ref_start,
                       ref_end => $ref_end,

                       evalue => $evalue,
                       bitscore => $bitscore,
                   };

        push (@alignments, $struct);
    }

    return(@alignments);

}




####
sub get_seq_lengths {
    my ($refseqs_fa) = @_;

    my %lengths;
    my $acc;
    
    open (my $fh, $refseqs_fa) or die $!;
    while (<$fh>) {
        if (/>(\S+)/) {
            $acc = $1;
        }
        else {
            my $seq = $_;
            chomp $seq;
            my $len = length($seq);
            $lengths{$acc} += $len;
        }
    }
    close $fh;
    
    return(%lengths);
}
    
####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    return;
}
