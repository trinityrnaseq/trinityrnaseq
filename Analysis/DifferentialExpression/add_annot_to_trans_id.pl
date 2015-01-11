#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 annot.gff3 target_file.txt\n\n";

my $annot_gff3 = $ARGV[0] or die $usage;
my $target_file = $ARGV[1] or die $usage;


main: {

    my %trans_id_to_annot = &parse_trans_names($annot_gff3);

    open (my $fh, $target_file) or die "Error, cannot open file $target_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        if (my $name = $trans_id_to_annot{ $x[0] }) {
            $x[0] .= " $name";
            $x[0] =~ s/\s/_/g;
        }

        print join("\t", @x) . "\n";
    }
    close $fh;



    exit(0);
    
}


####
sub parse_trans_names {
    my ($gff3_file) = @_;


    my %gene_to_name;
    my %mRNA_to_gene;

    open (my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
    while (<$fh>) {
        unless (/\w/) { next; }
        
        my @x = split(/\t/);
        my $feat_type = $x[2];
        my $info = $x[8];
        
        if ($info =~ /ID=([^;]+)/) {
            my $id = $1;

            if ($feat_type eq "mRNA") {
                $info =~ /Parent=([^;]+)/ or die "Error, cannot get parent info from $_";
                my $gene = $1;
                $mRNA_to_gene{$id} = $gene;
            }
            elsif ($feat_type eq "gene") {
                $info =~ /Name=(.*)/ or die "Error, cannot get name from $_";
                my $name = $1;
                $gene_to_name{$id} = $name;
            }
        }
    }
    close $fh;

    my %trans_to_name;
    foreach my $mRNA (keys %mRNA_to_gene) {
        my $gene = $mRNA_to_gene{$mRNA};
        my $name = $gene_to_name{$gene};
        $trans_to_name{$mRNA} = $name;
    }


    return(%trans_to_name);
}


        
