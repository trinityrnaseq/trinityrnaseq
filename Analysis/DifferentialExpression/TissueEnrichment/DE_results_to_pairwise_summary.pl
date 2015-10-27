#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $usage = "usage: $0 sample_avg_expr.matrix edgeR_directory/ [FDR=0.05]\n\n";

my $sample_expr_matrix = $ARGV[0] or die $usage;
my $DE_dir = $ARGV[1] or die $usage;
my $MAX_FDR = $ARGV[2];
unless (defined $MAX_FDR) {
    $MAX_FDR = 0.05;
}

my @DE_result_files = <$DE_dir/*.DE_results>;
unless (@DE_result_files) {
    die "Error, cannot find \*.DE_results files at $DE_dir ";
}

main: {
    

    my %gene_to_sample_expr_val = &parse_expression_matrix($sample_expr_matrix);

    print join("\t", "#feature", "sample_A", "sample_B", "log2(exprA)", "log2(exprB)", "logFC", "FDR") . "\n";
    
    foreach my $DE_result_file (@DE_result_files) {
        print STDERR "-processing DE file: $DE_result_file\n";
        
        $DE_result_file =~ /\.([^\.\/]+)_vs_([^\.\/]+).[^\.]+.DE_results/ or die "Error, cannot parse filename: $DE_result_file";
        my $sample_A = $1;
        my $sample_B = $2;
        
        open (my $fh, $DE_result_file) or die $!;
        my $header = <$fh>;
        chomp $header;
        my @header_fields = split(/\t/, $header);
        my $FDR_field = undef;

        my $line_counter = 0;
        while(<$fh>) {
            $line_counter++;
            chomp;
            my @x = split(/\t/);

            # match up header with data fields.
            if ($line_counter == 1) {
                if (scalar(@header_fields) == scalar(@x) -1) {
                    unshift(@header_fields, 'id');
                }
                elsif (scalar(@header_fields) != scalar(@x)) {
                    die "Error, disconnect between header line and data line, number of fields are unequal and header isn't one short:\n"
                        . "header: @header_fields\n"
                        . "line: @x\n";
                }
                for (my $i = 0; $i <= $#header_fields; $i++) {
                    if ($header_fields[$i] eq 'FDR') {
                        $FDR_field = $i;
                    }
                }
                unless ($FDR_field) {
                    die "Error, couldn't identify column corresponding to FDR";
                }
            }
            
            my $feature = $x[0];
            my $feature_expr = $gene_to_sample_expr_val{$feature};
            unless (defined $feature_expr) {
                die "Error, no expression values stored for [$feature] ";
            }
            

            my $FDR = $x[$FDR_field];
            if ($FDR <= $MAX_FDR) {
                my $expr_sample_A = $feature_expr->{$sample_A};
                my $expr_sample_B = $feature_expr->{$sample_B};
                
                unless (defined $expr_sample_A && defined $expr_sample_B) {
                    die "Error, no expr value for feature: $feature, $sample_A [$expr_sample_A] or $sample_B [$expr_sample_B] " . Dumper($gene_to_sample_expr_val{$feature});
                }
                my $log_expr_sample_A = log($expr_sample_A+1)/log(2);
                my $log_expr_sample_B = log($expr_sample_B+1)/log(2);
                
                my $log_FC = sprintf("%.2f", $log_expr_sample_A - $log_expr_sample_B);

                print join("\t", $feature, $sample_A, $sample_B, $log_expr_sample_A, $log_expr_sample_B, $log_FC, $FDR) . "\n";
            }
        }
    }

    print STDERR "\nDone\n\n";

    exit(0);

}


####
sub parse_expression_matrix {
    my ($expr_matrix_file) = @_;

    print STDERR "\nReading matrix: $expr_matrix_file ... ";

    my $cmd = "wc -l $expr_matrix_file ";
    
    my $num_lines = `$cmd`;
    if ($?) {
        die "Error, cmd; $cmd died with ret $?";
    }
    $num_lines =~ /(\d+)/;
    $num_lines = $1 or die "Error, cannot count number of lines from: $num_lines, cmd: $cmd";
    
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

        #unless ($feature_name eq "c1088792_g2_i5^sp|Q06441|TSP4_XENLA^COMP^sigP") { next; }
        
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
        #if ($counter > 10) { last; } # debug
    }

    close $fh;

    print STDERR "\n\nDone reading matrix.\n";

    #print Dumper(\%gene_to_sample_expr_val);
    
    return(%gene_to_sample_expr_val);
}

