#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename;

my $usage = "usage: $0 fileA.results.txt fileB.results.txt\n\n";


my @files = @ARGV;
unless (@files) {
    die $usage;
}



my %data;


foreach my $file (@files) {
    
    open (my $fh, $file) or die $!;
    $file = basename($file);
    my $header = <$fh>;
    while (<$fh>) {
        chomp;
        my ($acc, $conc, $foldchange, $pval, $fdr) = split(/\t/);
        
        $data{$acc}->{$file} = { pval => sprintf("%.2e", $pval),
                                 fdr => sprintf("%.2e", $fdr), 
                                 fc => sprintf("%.1f", $foldchange),
                             };
    }
    
}

print "#acc";
foreach my $file (@files) {
    print "\t" . basename($file);
}
print "\tmost_signif_P\tcorresp_FDR\n";



foreach my $acc (keys %data) {

    print "$acc";

    my $most_signif_entry;
    
    foreach my $file (@files) {
        my $basefile = basename($file);
        
        if (my $struct = $data{$acc}->{$basefile}) {
            my $pval = $struct->{pval};
            my $fdr = $struct->{fdr};
            my $fc = $struct->{fc};
            print "\t(FC:$fc,P:$pval,FDR:$fdr)";
            if ( (! $most_signif_entry) || $most_signif_entry->{pval} > $pval) {
                $most_signif_entry = $struct;
            }
        }
        else {
            print "\tNA";
        }
    }
    
    print "\t" . $most_signif_entry->{pval} . "\t" . $most_signif_entry->{fdr} . "\n";
}


exit(0);

