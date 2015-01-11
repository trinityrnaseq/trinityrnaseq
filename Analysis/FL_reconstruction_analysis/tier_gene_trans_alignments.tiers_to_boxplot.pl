#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

my $usage = "usage: $0 fileA.tiers [fileB.tiers ...]\n\n";

my @tier_files = @ARGV or die $usage;



main: {

    foreach my $tier_file (@tier_files) {
        
        print STDERR "// processing $tier_file\n";
        
        my @top_tier_per_gene = &get_top_tier_per_gene($tier_file);
        
        open (my $ofh, ">$tier_file.top_tiers") or die $!;
        print $ofh join("\n", @top_tier_per_gene) . "\n";
        close $ofh;
        
    }

    ## generate boxplot
    open (my $ofh, ">__tmp_tiers_boxplot.R") or die $!;
    my @name_tokens;
    foreach my $tier_file (@tier_files) {
        my $token = basename($tier_file);
        $token =~ s/\W/_/g;
        push (@name_tokens, $token);
        print $ofh "$token = read.table(\"$tier_file.top_tiers\", header=F);\n";
    }
    print $ofh "pdf(file=\"tiers.boxplot.pdf\");\n";
    print $ofh "boxplot(c(" . join(",", @name_tokens) . "), names=c(\"" . join("\",\"", @name_tokens) . "\"), outline=T, las=2);\n";
    print $ofh "dev.off();\n";
    close $ofh;
    
    system("R --vanilla -q < __tmp_tiers_boxplot.R");

    exit($?);
}

####
sub get_top_tier_per_gene {
    my ($tier_file) = @_;

    my %gene_to_top_tier;

    open (my $fh, $tier_file) or die $!;
    while (<$fh>) {
        #print;
        
        chomp;
        
        if (/^Tier\[(\d+)\]\s+(\S+)/) {
            my $tier_val = $1;
            my $gene_id = $2;
            if ( (! exists $gene_to_top_tier{$gene_id}) 
                 ||
                 $gene_to_top_tier{$gene_id} < $tier_val) {
                
                $gene_to_top_tier{$gene_id} = $tier_val;
                
                #print STDERR "$gene_id => $tier_val\t" . substr($_, 0, 100) . "\n";
            }
        }
        else {
            die "Error, cannot parse $_";
        }
        
    }
    close $fh;

    my @vals = values %gene_to_top_tier;

    return(@vals);
}
