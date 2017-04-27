#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Process_cmd;


my $usage = "\n\n\tusage: $0 Trinity.fasta [out_prefix='GC_content']\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $out_prefix = $ARGV[1] || "GC_content";

main: {

    my $dat_out_file = "$out_prefix.dat";

    open(my $ofh, ">$dat_out_file") or die "Error, cannot write to $dat_out_file";
    
    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();

        my $seq_len = length($sequence);

        my $gc_count = 0;
        while ($sequence =~ /[gc]/ig) {
            $gc_count++;
        }
        my $pct_gc = sprintf("%.2f", $gc_count / $seq_len * 100);

        print $ofh join("\t" ,$acc, $pct_gc) . "\n";
    }

    # generate histogram
    my $R_code = <<__eoR__;

    data = read.table("$dat_out_file", header=F, row.names=1)
    pdf("$dat_out_file.hist.pdf")
    hist(data[,1], br=100)
    message("\n\nmean: ", sprintf("%.2f", mean(data[,1])), ", median: ", median(data[,1]), "\n\n")
    dev.off()

__eoR__
    
;
    
    {
        my $Rscript_file = "$out_prefix.R";
        open (my $ofh, ">$Rscript_file") or die "Error, cannot write to file: $Rscript_file";
        print $ofh $R_code;
        close $ofh;

        # run it
        my $cmd = "Rscript $Rscript_file";
        &process_cmd($cmd);
    }
    
    exit(0);
}
