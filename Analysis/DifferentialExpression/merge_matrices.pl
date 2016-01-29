#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin;
use File::Basename;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Data::Dumper;


my @matrices = @ARGV;

unless (scalar @matrices > 1) {
    die "\n\n\tusage: $0 matrixA matrixB ...\n\n";
}

my %matrix;
my %genes;

main: {
    
    foreach my $matrix (@matrices) {

        &parse_matrix($matrix);
        
    }

    ## output new matrix:

    my @colnames = sort keys %matrix;
    print "\t" . join("\t", @colnames) . "\n";

    foreach my $gene (sort keys %genes) {
        
        print "$gene";
        foreach my $colname (@colnames) {
            my $val = $matrix{$colname}->{$gene};
            unless (defined $val) {
                $val = "NA";
            }
            print "\t$val";
        }
        print "\n";
    }
    

    exit(0);

}




####
sub parse_matrix {
    my ($matrix_file) = @_;
    
    open (my $fh, $matrix_file);
    my $header = <$fh>;
    chomp $header;
    my @pos_to_col = split(/\t/, $header);
    my $check_column_ordering_flag = 0;
    
    foreach my $sample (@pos_to_col) {
        if (exists $matrix{$sample}) {
            die "Error, already encountered column header: $sample, cannot have redundant column names across matrices.";
        }
    }

    
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        
        unless ($check_column_ordering_flag) {
            if (scalar(@x) == scalar(@pos_to_col) + 1) {
                ## header is offset, as is acceptable by R
                ## not acceptable here.  fix it:
                unshift (@pos_to_col, "");
            }
            $check_column_ordering_flag = 1;
          
        }
        

        my $gene = $x[0];
        $genes{$gene} = 1;
        
        for (my $i = 1; $i <= $#x; $i++) {
            my $col = $pos_to_col[$i];
            my $val = $x[$i];
                
            $matrix{$col}->{$gene} = $val;
            
        }
        
    }
   

    return %matrix;

}

