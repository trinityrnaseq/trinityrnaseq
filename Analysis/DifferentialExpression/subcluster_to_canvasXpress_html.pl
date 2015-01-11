#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::Bin/../TrinityWeb/cgi-bin/PerlLib/");
use CanvasXpress::Heatmap;
use CanvasXpress::PlotOnLoader;
use CanvasXpress::Line;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use CGI;


my $usage = <<__EOUSAGE__;

#############################################################################################
#
# Required:
#
#  --matrix <string>     matrix (best if log2, centered)
#
# Optional:
#
#  --igv_coords <string>   text file, formatted: feature_name (tab) chr (tab) lend (tab) rend
#
#  --cluster_samples       cluster the samples in the heatmap
#
##############################################################################################


__EOUSAGE__

    ;


my $matrix;
my $igv_coords_file;
my $sample_clusters_flag = 0;
my $help_flag;

&GetOptions('matrix=s' => \$matrix,
            'igv_coords=s' => \$igv_coords_file,
            'cluster_samples' => \$sample_clusters_flag,
            'help|h' => \$help_flag,
    );

if ($help_flag) {
    die $usage;
}

unless ($matrix) {
    die $usage;
}



main: {
    
    my $cgi = new CGI();
    #print $cgi->header();

    my $plot_loader_func_name = "plot_loader_$$";
    
    print $cgi->start_html(-title => 'html matrix',
                           -onLoad => $plot_loader_func_name . "();",
     );

    my $plot_loader = new CanvasXpress::PlotOnLoader($plot_loader_func_name);
    
    open (my $fh, $matrix) or die $!;
    my $header = <$fh>;
    chomp $header;
    $header =~ s/^\s+//;

    my %accs;

    my @samples = split(/\s+/, $header);
    my @values;
    while (<$fh>) {
        chomp;
        my @x = split(/\s+/);
        push (@values, [@x]);
       
        my $acc = $x[0];
        $accs{$acc} = 1;
    }
    close $fh;

    my %inputs = ( samples => \@samples,
                   value_matrix => \@values,
                   cluster_features => 1,
                   cluster_samples => $sample_clusters_flag,
        );
    
    if ($igv_coords_file) {
        &add_IGV_links($igv_coords_file, \%accs);
        
        &add_click_event(\%inputs);

    }
    

    my $heatmap_obj = new CanvasXpress::Heatmap("heatmap_$$");
    print $heatmap_obj->draw(%inputs, 'dendrogramSpace' => 1);
    $plot_loader->add_plot($heatmap_obj);
    
    my $line_graph_obj = new CanvasXpress::Line("line_$$");
    print $line_graph_obj->draw(%inputs, 'graphOrientation' => 'vertical');
    $plot_loader->add_plot($line_graph_obj);
    
    print $plot_loader->write_plot_loader();
        
    print $cgi->end_html();
    
    exit(0);
}

####
sub add_IGV_links {
    my ($coords_file, $accs_href) = @_;

    print "<script>\n"
        . "   var igv_lookup = {};\n";

    open (my $fh, $coords_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($acc, $chr, $lend, $rend) = split(/\t/);

        ## offset by 10%
        my $delta = $rend - $lend + 1;
        my $offset = int(0.1 * $delta);
        $lend -= $offset;
        if ($lend < 1) {
            $lend = 1;
        }
        $rend += $offset;
        
        unless ($accs_href->{$acc}) { next; }
        if ($acc =~ /\'/) { next; } # sorry, aint gonna work
        my $val = "$chr:$lend-$rend";
                
        print "  igv_lookup[\'$acc\'] = \'$val\';\n";
    }
    close $fh;
    
    print "</script>\n";

    return;
}

####
sub add_click_event {
    my ($inputs_href) = @_;

    $inputs_href->{events} = { 'click' => "var gene = o['y']['vars'][0];\n"
                                   . "document.location.href=\'http://localhost:60151/goto?locus=\' + igv_lookup[gene];\n" };

    return;
}
