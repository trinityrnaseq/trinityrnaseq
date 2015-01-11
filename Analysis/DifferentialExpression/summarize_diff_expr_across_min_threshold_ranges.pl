#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 all_diff_expression_results.txt\n\n";

my $results_summary_file = $ARGV[0] or die $usage;

my @signif_results;
{
    open (my $fh, $results_summary_file) or die "Error, cannot open file $results_summary_file";
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my @x = split(/\t/);
        
        my ($sampleA, $sampleB, $acc, $logFC, $logCPM, $pvalue, $fdr) = @x;

        my $struct = { sampleA => $sampleA,
                       sampleB => $sampleB,
                       acc => $acc,
                       logFC => $logFC,
                       pvalue => $pvalue,
                       fdr => $fdr,
        };

        push (@signif_results, $struct);
    }
    close $fh;

}


print "<html>\n";
print "<table>\n";

print "<tr><th>.</th><th>" . join("</th><th>1e-", (1..10)) . "</th></tr>\n";


for my $min_logfc (1..10) {

    
    my @counts;
    
    print "<tr><th>" . 2**$min_logfc . "</th>";

    for my $min_fdr (1..10) {

        my %accs = map { $_->{acc} => 1 }  grep { $_->{fdr} <= 10**(-1*$min_fdr) && $_->{logFC} >= $min_logfc } @signif_results;

        my $count = scalar(keys %accs);
        #print join("\t", $min_logfc, $min_fdr, $count) . "\n";
        push (@counts, $count);
    }

    print "<td>" . join("</td><td>", @counts) . "</td><tr>\n";
    
}
print "</table>\n";
print "</html>\n";

exit(0);

