#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#####################################################
#
#  --DE_results <string>          DE_results file
#
#  --gene_trans_map <string>      gene-trans-map file
#
#####################################################


__EOUSAGE__

    ;

my $help_flag;


my $DE_results_file;
my $gene_trans_map_file;

&GetOptions ( 'h' => \$help_flag,

              'DE_results=s' => \$DE_results_file,
              'gene_trans_map=s' => \$gene_trans_map_file,
              
    );


if ($help_flag) {
    die $usage;
}

unless ($DE_results_file && $gene_trans_map_file) {
    die $usage;
}


main: {

    my %trans_to_gene = &parse_gene_trans_map($gene_trans_map_file);

    my %gene_to_up_down_trans;
    
    open(my $fh, $DE_results_file) or die "Error, cannot open file $DE_results_file";
    my $header = <$fh>;
    chomp $header;
    my @header_fields = split(/\t/, $header);
    my %header_col_index;
    
    while (<$fh>) {
        chomp;
        my @fields = split(/\t/);
        unless (%header_col_index) {
            if (scalar(@fields) == scalar(@header_fields) + 1) {
                unshift(@header_fields, "");
            }
        }
        
        my $trans_name = $fields[0];

        my $gene_name = $trans_to_gene{$trans_name} or die "Error, cannot find gene_id corresponding to trans [$trans_name] ";
        
        my $struct = { name => $trans_name };
        for (my $i = 1; $i <= $#fields; $i++) {
            my $col_header = $header_fields[$i];
            my $val = $fields[$i];
            $struct->{$col_header} = $val;
        }
        my $fdr = $struct->{FDR};
        unless (defined $fdr) { 
            die "Error, no FDR set for @fields";
        }
        
        unless ($fdr <= 0.05) { next; }

        my $logFC = $struct->{logFC};
        unless (defined $logFC) {
            die "Error, no logFC set for @fields";
        }
        unless (abs($logFC) >= 1) { # 2-fold min
            next; 
        }

        if ($logFC > 0) {
            push (@{$gene_to_up_down_trans{$gene_name}->{UP}}, $struct);
        }
        else {
            push(@{$gene_to_up_down_trans{$gene_name}->{DOWN}}, $struct);
        }
        
        
    
    }
    close $fh;

    
    foreach my $gene_name (keys %gene_to_up_down_trans) {
        my @up;
        if (exists $gene_to_up_down_trans{$gene_name}->{UP}) {
            @up = @{$gene_to_up_down_trans{$gene_name}->{UP}};
        }
        
        my @down;
        if (exists $gene_to_up_down_trans{$gene_name}->{DOWN}) {
            @down = @{$gene_to_up_down_trans{$gene_name}->{DOWN}};
        }
        
        if (@up && @down) {
            ## got DTU candidate:
            print "# $gene_name\n";
            foreach my $up_entry (@up) {
                print join("\t", "UP", &struct_to_string($up_entry)) . "\n";
            }
            foreach my $down_entry (@down) {
                print join("\t", "DOWN", &struct_to_string($down_entry)) . "\n";
            }
            print "\n"; # spacer
        }
    }

    exit(0);
    

}


####
sub parse_gene_trans_map {
    my ($gene_trans_map_file) = @_;

    my %trans_to_gene;

    open(my $fh, $gene_trans_map_file) or die "Error, cannot open file $gene_trans_map_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
            
        my ($gene, $trans) = split(/\t/);
        
        $trans_to_gene{$trans} = $gene;
    }

    return(%trans_to_gene);
}
        
####
sub struct_to_string {
    my ($struct) = @_;

    my $ret_text = join("\t", $struct->{name}, $struct->{sampleA}, $struct->{sampleB}, 
                        $struct->{logFC}, $struct->{FDR});
    
    return($ret_text);
}

