#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use lib "$FindBin::RealBin/../../PerlLib";
use Pipeliner;
use File::Basename;


my $usage = <<__EOUSAGE__;

#######################################################################################
#
#  --outfmt6_grouped <string>      outfmt6 grouped output
#
#  --min_pct_len <int>             minimum percent length covered by pairwise matches
#
#  --min_per_id <int>              minimum percent identity
#
#  --inflation_factor <float>      inflation factor for MCL clustering
#
#######################################################################################


__EOUSAGE__

    ;



my $help_flag;
my $outfmt6_grouped_file;
my $inflation_factor;
my $min_pct_len;
my $min_per_id;

&GetOptions ( 'h' => \$help_flag,
              
              'outfmt6_grouped=s' => \$outfmt6_grouped_file,

              'min_pct_len=i' => \$min_pct_len,

              'min_per_id=i' => \$min_per_id,
              
              'inflation_factor=f' => \$inflation_factor,
              
    );

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, dont understand parameters @ARGV";
}


unless ($outfmt6_grouped_file && $inflation_factor && $min_pct_len && $min_per_id) {
    die $usage;
}


# add MCL to PATH setting
$ENV{PATH} = "/seq/regev_genome_portal/SOFTWARE/MCL/bin/:$ENV{PATH}";



main: {

    
    my $filtered_hits = basename($outfmt6_grouped_file . ".minLEN_${min_pct_len}_pct_len.minPID_${min_per_id}.abc");
    my $checkpoint = ".$filtered_hits.ok";
    if (! -e $checkpoint) {
        
        my %best_hits;
        


        open (my $fh, $outfmt6_grouped_file) or die "Error, cannot open file $outfmt6_grouped_file";
        while (<$fh>) {
            if (/^\#/) { next; }
            chomp;
            my @x = split(/\t/);
            my ($transA, $transB, $per_id, $E_value, @rest) = split(/\t/);
            my $per_len_match = pop @rest;
            
            if ($per_len_match >= $min_pct_len && $per_id >= $min_per_id) {
                my $geneA = &parse_gene_name($transA);
                my $geneB = &parse_gene_name($transB);
                
                if ($geneA eq $geneB) { next; } 
                
                ($geneA, $geneB) = sort ($geneA, $geneB);
                my $gene_pair_token = join("$;", $geneA, $geneB);
               

                my $lowest_evalue = $best_hits{$gene_pair_token};
                if ( (! defined $lowest_evalue) || $lowest_evalue > $E_value) {
                    $best_hits{$gene_pair_token} = $E_value;
                }
                
                
            }
        }

        # write best hits file
        open (my $ofh, ">$filtered_hits") or die "Error, cannot write to $filtered_hits";        
        foreach my $gene_pair_token (keys %best_hits) {
            my ($geneA, $geneB) = split(/$;/, $gene_pair_token);
            my $E_value = $best_hits{$gene_pair_token};
            print $ofh join("\t", $geneA, $geneB, $E_value) . "\n";
        }
        close $ofh;
        
        `touch $checkpoint`;
    }

    my $pipeliner = new Pipeliner(-verbose => 1);
    my $cmd = "mcxload -abc $filtered_hits --stream-mirror --stream-neg-log10 "
        . " -stream-tf 'ceil(200)' -o $filtered_hits.mci -write-tab $filtered_hits.tab";
    
    $pipeliner->add_commands( new Command($cmd, ".$filtered_hits.tab.ok") );
    
    $inflation_factor = sprintf("%.1f", $inflation_factor);
    my $inflation_factor_dec_removed = $inflation_factor;
    $inflation_factor_dec_removed =~ s/\.//;

    $cmd = "mcl $filtered_hits.mci -I $inflation_factor";
    my $mcl_outfile = "out.$filtered_hits.mci.I$inflation_factor_dec_removed";
    $pipeliner->add_commands( new Command($cmd, ".$mcl_outfile.ok"));
    
    
    $cmd = "mcxdump -icl $mcl_outfile -tabr $filtered_hits.tab -o dump.$mcl_outfile";
    $pipeliner->add_commands( new Command($cmd, ".dump.$mcl_outfile.ok"));


    $pipeliner->run();


    exit(0);

}

####
sub parse_gene_name {
    my ($trans_info) = @_;

    my ($gene_symbol, $trans_id);
    if (/;/) {
        ($trans_id, $gene_symbol) = split(/;/, $trans_info);
    }
    elsif (/\|/) {
        ($gene_symbol, $trans_id) = split(/\|/, $trans_info);
    }
    

    if ($gene_symbol) {
        return($gene_symbol);
    }
    else {
        return($trans_info);
    }
}
