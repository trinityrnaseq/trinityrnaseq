#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 sample_avg_expr.matrix EBSeq_directory/ [minPosteriorProb=0.95]\n\n";

my $sample_expr_matrix = $ARGV[0] or die $usage;
my $EBSeq_dir = $ARGV[1] or die $usage;
my $min_posterior = $ARGV[2];
if (! defined $min_posterior) {
    $min_posterior = 0.95;
}

my @DE_result_files = <$EBSeq_dir/*.EBSeq>;
unless (@DE_result_files) {
    die "Error, cannot find \*.EBSeq results files at $EBSeq_dir ";
}

main: {
    

    my %gene_to_sample_expr_val = &parse_expression_matrix($sample_expr_matrix);

    foreach my $DE_result_file (@DE_result_files) {

        $DE_result_file =~ /\W([^\.\/]+)-vs-([^\.\/]+).\w+\.EBSeq/ or die "Error, cannot parse filename: $DE_result_file";
        my $sample_A = $1;
        my $sample_B = $2;

        open (my $fh, $DE_result_file) or die $!;
        my $header = <$fh>;
        while(<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $feature = $x[0];
            $feature =~ s/\"//g;
            my $posterior = $x[2];
            if ($posterior >= $min_posterior) {
                my $expr_sample_A = $gene_to_sample_expr_val{$feature}->{$sample_A};
                my $expr_sample_B = $gene_to_sample_expr_val{$feature}->{$sample_B};
                
                unless (defined $expr_sample_A && defined $expr_sample_B) {
                    die "Error, no expr value for feature: $feature, $sample_A [$expr_sample_A] or $sample_B [$expr_sample_B] " . Dumper($gene_to_sample_expr_val{$feature});
                }
                my $log_expr_sample_A = log($expr_sample_A+1)/log(2);
                my $log_expr_sample_B = log($expr_sample_B+1)/log(2);
                
                my $log_FC = sprintf("%.2f", $log_expr_sample_A - $log_expr_sample_B);

                print join("\t", $feature, $sample_A, $sample_B, $log_expr_sample_A, $log_expr_sample_B, $log_FC, $posterior) . "\n";
            }
        }
    }


    exit(0);

}


####
sub parse_expression_matrix {
    my ($expr_matrix_file) = @_;

    print STDERR "\nReading matrix: $expr_matrix_file ... ";

    my $num_lines = `wc -l $expr_matrix_file | cut -f1 -d ' '`;
    chomp $num_lines;
    
    print STDERR " $num_lines rows of matrix detected.\n\n";
    
    my %gene_to_sample_expr_val;

    open (my $fh, $expr_matrix_file) or die "Error, cannot open file $expr_matrix_file";
    my $header = <$fh>;
    chomp $header;
    $header =~ s/^\s+//;
    my @sample_names = split(/\t/, $header);
    
    my $counter = 0;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $feature_name = shift @x;
        
        unless (scalar @x == scalar @sample_names) {
            die "Error, number of samples: " . scalar (@sample_names) . " doesn't match number of values read: " . scalar(@x) . "  ";
        }

        for (my $i = 0; $i <= $#sample_names; $i++) {
            my $sample = $sample_names[$i];
            my $val = $x[$i];
            
            $gene_to_sample_expr_val{$feature_name}->{$sample} = $val;
        }
        

        $counter++;
        if ($counter % 10000 == 0) {
            my $pct_done = sprintf("%.2f", $counter/$num_lines * 100);
            print STDERR "\r[$pct_done %] matrix read.   ";
        }
    }

    close $fh;

    print STDERR "\n\nDone reading matrix.\n";

    return(%gene_to_sample_expr_val);
}

