#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use FindBin;
use File::Basename;

use lib ("$FindBin::Bin/PerlLib");

use CanvasXpress::Heatmap;
use BHStats;

$|++;

my $cgi = new CGI();

print $cgi->header();

my %params = $cgi->Vars();

print $cgi->start_html('-title' => 'Diff Expression Analysis Using CanvasXpress');

my $data_dir = $params{'data_dir'};

## defaults for params
{
    unless ($data_dir) {
        $data_dir = $params{'data_dir'} = 'sample_data';
        print "<p><b>Sample data shown</b></p>\n";
    }
    
    unless (exists $params{max_FDR}) {
        $params{max_FDR} = 0.001;
    }
    unless (exists $params{min_FC}) {
        $params{min_FC} = 4;
    }
    
}



my $max_FDR = $params{max_FDR};
my $min_FC = $params{min_FC};

my $FPKM_matrix_file = "$data_dir/matrix.TMM_normalized.FPKM";
my $diff_expressed_results_file = "$data_dir/all_diff_expression_results.txt";

unless (-s $FPKM_matrix_file) {
    die "Error, cannot find file: $FPKM_matrix_file";
}
unless (-s $diff_expressed_results_file) {
    die "Error, cannot find file $diff_expressed_results_file";
}

my $MAX_SHOW = 200;


main: {

    my %diff_expressed_transcripts = &get_diff_expressed_transcripts($diff_expressed_results_file, $min_FC, $max_FDR);

    my $count_diff_expr = scalar(keys %diff_expressed_transcripts);

    if ($count_diff_expr > $MAX_SHOW) {
        my @accs = reverse sort {$diff_expressed_transcripts{$a} <=> $diff_expressed_transcripts{$b}} keys %diff_expressed_transcripts;
        @accs = @accs[0..$MAX_SHOW-1];
        my %adj_diff_expr_trans = map { $_ => $diff_expressed_transcripts{$_} } @accs;
        %diff_expressed_transcripts = %adj_diff_expr_trans;
    }
    
    my @matrix_entries = &get_matrix_entries($FPKM_matrix_file, \%diff_expressed_transcripts);

    my $samples_aref = shift @matrix_entries;

    my @log2_median_centered_matrix = &log2_median_center(@matrix_entries);

    

    
        
    

    my $heatmap_html = &CanvasXpress::Heatmap::draw( samples => $samples_aref,
                                                     value_matrix => \@log2_median_centered_matrix,
                                                     
                                                     cluster_features => $params{cluster_transcripts},
                                                     cluster_samples => $params{cluster_samples},
                                                     
                                                     );
    

    
    my $controls_html = &write_heatmap_controls_html(\%params);

    print $controls_html;
    
    print "<p><b>$count_diff_expr</b> transcripts identified as differentially expressed given thresholds.</p>\n";
    if ($count_diff_expr > $MAX_SHOW) {
        print "<p>Only the top $MAX_SHOW most statistically significant transcripts are shown.</p>\n";
    }
    
    print $heatmap_html;
    

}


print $cgi->end_html();


exit(0);

####
sub get_diff_expressed_transcripts {
    my ($diff_expressed_results_file, $min_FC, $max_FDR) = @_;

    my $min_logFC = log($min_FC)/log(2);

    my %diff_trans;

    open (my $fh, $diff_expressed_results_file) or die $!;
    while (<$fh>) {
        chomp;
        if (/^\#/) { next; }
        my ($sampleA, $sampleB, $transcript, $logConc, $logFC, $Pvalue, $FDR) = split(/\t/);
        
        if (abs($logFC) >= $min_logFC && $FDR <= $max_FDR) {
            
            if (exists $diff_trans{$transcript}) {
                if ($FDR < $diff_trans{$transcript}) {
                    $diff_trans{$transcript} = $FDR;
                }
            }
            else {
                $diff_trans{$transcript} = $FDR;
            }
        }
    }
    close $fh;
    

    return(%diff_trans);
}

####
sub get_matrix_entries {
    my ($matrix_file, $diff_expressed_trans_href) = @_;

    

    my @matrix_lines;
    
    open (my $fh, $matrix_file) or die "Error, cannot open file $matrix_file";
    my $header = <$fh>;
    chomp $header;
    $header =~ s/^\s+//;
    my @samples = split(/\t/, $header);
    
    push (@matrix_lines, [@samples]);

    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        if (exists $diff_expressed_trans_href->{$x[0]}) {
            push (@matrix_lines, [@x]);
        }
    }
    close $fh;


    return(@matrix_lines);
}

####
sub log2_median_center {
    my (@matrix_entries) = @_;
    
    my @adj_entries;

    foreach my $entry (@matrix_entries) {
        my ($feature_name, @vals) = @$entry;
        
        foreach my $v (@vals) {
            $v = log($v+1)/log(2);
        }
        my $median = &BHStats::median(@vals);
        foreach my $v (@vals) {
            $v -= $median;
        }
        
        push (@adj_entries, [$feature_name, @vals]);
    }

    return(@adj_entries);
}

        
####
sub write_heatmap_controls_html {
    my ($params_href) = @_;

    my $html = "<div id='control_box'>\n";
    
    $html .= "<form id='control_form' method='get' action='diff_express.cgi' >\n";
    $html .= "<ul>\n";
    $html .= "   <li>min_FC: <input type='text' name='min_FC' value=\'$params_href->{min_FC}\' />\n";
    $html .= "   <li>max_FDR: <input type='text' name='max_FDR' value=\'$params_href->{max_FDR}\' />\n";
    {
        $html .= "   <li><input type='checkbox' name='cluster_transcripts' value='1' ";
        if (exists $params_href->{cluster_transcripts}) {
            $html .= " checked ";
        }
        $html .= " />Cluster transcripts\n";
    }
    
    {
        $html .= "   <li><input type='checkbox' name='cluster_samples' value='1' ";
        if (exists $params_href->{cluster_samples}) {
            $html .= " checked ";
        }
        $html .= " />Cluster samples\n"; 
    }

    $html .= "   <li><input type='submit' />\n";
    $html .= "</ul>\n";

    $html .= "<input type='hidden' name='data_dir' value=\'$params_href->{data_dir}\' />\n";
    
    $html .= "</form>\n";
    $html .= "</div>\n";


    return($html);
}

