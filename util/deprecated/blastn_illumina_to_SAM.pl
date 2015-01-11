#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use SAM_entry;
use Nuc_translator;

my $usage = "usage: $0 reads.fasta target.fasta\n\n";

my $reads_fasta = $ARGV[0] or die $usage;
my $target_fasta = $ARGV[1] or die $usage;

main: {
    
    my $cmd = "fasta_to_tab.pl $reads_fasta > $reads_fasta.tab";
    &process_cmd($cmd) unless (-s "$reads_fasta.tab");

    $cmd = "sort -k1,1 $reads_fasta.tab > $reads_fasta.tab.sort";
    &process_cmd($cmd) unless (-s "$reads_fasta.tab.sort");

    $cmd = "formatdb -i $target_fasta -p F";
    &process_cmd($cmd) unless (-s "$target_fasta.nsq");
    
    $cmd = "blastall -p blastn -d $target_fasta -i $reads_fasta -G -1 -E -1 -r 4 -q -4 -m 8 -v 1 -b 1 -e 1 > blastn.m8";
    &process_cmd($cmd) unless (-s "blastn.m8");

    $cmd = "sort -k 1,1 blastn.m8 > blastn.m8.sort";
    &process_cmd($cmd) unless (-s "blastn.m8.sort");

    
    my $sam_file = &convert_to_sam("blastn.m8.sort", "$reads_fasta.tab.sort");

    ## convert to bam
    $cmd = "samtools faidx $target_fasta";
    &process_cmd($cmd) unless (-s "$target_fasta.fai");
    

    my $bam_file = $sam_file;
    $bam_file =~ s/\.sam$//;
    
    $cmd = "samtools view -bt $target_fasta.fai $sam_file > $bam_file.bam";
    &process_cmd($cmd);
    
    $cmd = "samtools sort $bam_file.bam $bam_file";
    &process_cmd($cmd);

    $cmd = "samtools index $bam_file.bam";
    &process_cmd($cmd);
    
    
    exit(0);

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


####
sub convert_to_sam {
    my ($blast_sort, $reads_sort) = @_;

    my $sam_file = "$reads_sort.sam";
    open (my $sam_ofh, ">$sam_file") or die $!;

    my $prev = "";
   
    my $read_reader = Reads_reader->new($reads_sort);
    

    open (my $fh, $blast_sort) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $read_acc = $x[0];
        
        if ($read_acc eq $prev) { next; }
        $prev = $read_acc;
        
        my $read_seq = $read_reader->get_read($read_acc);

        

        my $target_acc = $x[1];
        my $lend = $x[6];
        my $rend = $x[7];
        my $m_lend = $x[8];
        my $m_rend = $x[9];
        
        my $read_length = length($read_seq);

        my $orient = '+';
        if ($m_rend < $m_lend) {
            $orient = '-';
            ($m_lend, $m_rend) = ($m_rend, $m_lend);
            
            # revcomp everything
            $read_seq = &reverse_complement($read_seq);
            $lend = $read_length - $lend + 1;
            $rend = $read_length - $rend + 1;
        
            ($lend, $rend) = ($rend, $lend);
        }
        
        my $match_len = $rend - $lend + 1;

        my $cigar_string = "${match_len}M";
        if ($lend != 1) {
            $cigar_string = ($lend-1) . "S" . $cigar_string;
        }
        if ($rend != $read_length) {
            my $delta = $read_length - $rend;
            $cigar_string .= $delta . "S";
        }
        
        my $flag = ($orient eq '-') ? 0x0010 : 0;
        
        my $sam_line = join("\t", $read_acc, $flag, $target_acc, $m_lend, 255, $cigar_string, '*', 0, 0, 
                            $read_seq, "#" x $read_length);
        
        
        print $sam_ofh $sam_line . "\n";
        

    }
    
    close $sam_ofh;
    

    return ($sam_file);
}





#####
package Reads_reader;

use strict;
use warnings;

sub new {
    my ($packagename) = shift;
    my $file = shift;

    open (my $fh, $file) or die $!;
    my $self = {

        fh => $fh,
        
    };

    bless($self, $packagename);

    return($self);
}
        

sub get_read {
    my $self = shift;
    my $acc = shift;
    
    my $fh = $self->{fh};

    while (<$fh>) {
        chomp;
        my ($read_acc, $seq_acc) = split(/\t/);
        
        if ($read_acc eq $acc) {
            return($seq_acc);
        }

        elsif ($read_acc gt $acc) {
            die "Error, must have missed acc: $acc in the sorted list.";
        }
    }
    

    die "Error, reached end of sorted read file without finding acc: $acc";
}



    

