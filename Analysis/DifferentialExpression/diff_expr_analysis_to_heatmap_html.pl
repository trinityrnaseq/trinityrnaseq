#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use FindBin;

use lib ("$FindBin::RealBin/../../PerlLib");
use CanvasXpress::Heatmap;


use CGI;
use CGI::Carp qw(fatalsToBrowser);


my $usage = <<__EOUSAGE__;

###################################################################################
#
# -R <string>  the filename for the stored RData (file.all.RData)
#
###################################################################################


__EOUSAGE__

    ;


my $help_flag = 0;
my $R_data_file;

&GetOptions ( 'h' => \$help_flag,
              'R=s' => \$R_data_file,
              );


if ($help_flag) {
    die $usage;
}

unless ($R_data_file) {
    die $usage;
}

####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd $cmd died with ret $ret";
    }

    return;
}



main: {
    
    unless (-s $R_data_file) {
        die "Error, cannot find pre-existing R-session data as file: $R_data_file";
    }
    

    my $matrix_data = "$R_data_file.ordered_gene_matrix";
    my $gene_tree = "$R_data_file.gene_nodist_tree";

    my $gene_tree_text;

    if (! (-s $matrix_data && -s $gene_tree)) {
        
        # generate required files.
        
        my $R_script = "__tmp_write_heatmap_html.R";
        
        open (my $ofh, ">$R_script") or die "Error, cannot write to file $R_script";
        
        print $ofh "source(\"$FindBin::RealBin/R/get_cluster_info.R\")\n";
        print $ofh "get_cluster_info(\"$R_data_file\")\n";
        close $ofh;
        
        &process_cmd("R --vanilla -q --slave < $R_script");

        
        open (my $fh, $matrix_data) or die "Error, cannot open file $matrix_data";
        my $header = <$fh>;
        chomp $header;
        my @gene_ids;
        my $counter = 0;
        while (<$fh>) {
            chomp;
            my ($gene_id, @expr_vals) = split(/\t/);
            push (@gene_ids, $gene_id);
        }
        close $fh;
        
        $gene_tree_text = `cat $gene_tree`;
        $gene_tree_text =~ s/;\s+$//g;
        
        $gene_tree_text = &convert_ids($gene_tree_text, \@gene_ids);
        
        open ($ofh, ">$gene_tree.ids_converted") or die $!;
        print $ofh $gene_tree_text;
        close $ofh;

    
    }
    else {
        $gene_tree_text = `cat $gene_tree.ids_converted`;
    }
    
    #############################################
    ## Generate HTML Heatmap  ###################
    #############################################
    
    print "<html>\n";


    
        
    my @sample_ids;
    my @feature_matrix;

    open (my $fh, $matrix_data) or die "Error, cannot open file $matrix_data";
    my $header = <$fh>;
    chomp $header;
    @sample_ids = split(/\t/, $header);
    my $counter = 0;
    while (<$fh>) {
        chomp;
        my ($gene_id, @expr_vals) = split(/\t/);
                
        push (@feature_matrix, [$gene_id, @expr_vals]);
        
        $counter++;
        # if ($counter > 10) { last; }
        
    }
    close $fh;
    
    my %heatmap_inputs = ( samples => [@sample_ids],
                           value_matrix => [@feature_matrix],
                           feature_tree => $gene_tree_text,
                           );

    print &CanvasXpress::Heatmap::draw(%heatmap_inputs);

    



    
#print "cx.clusterVariables();\n";
#print "cx.clusterSamples();\n";
    
    


    print "</html>\n";


    exit(0);
}

####
sub convert_ids {
    my ($gene_tree_text, $gene_ids_aref) = @_;

    my $counter = 0;
    foreach my $gene_id (@$gene_ids_aref) {
        $counter++;

        $gene_tree_text =~ s/([\(,])$counter([\),])/$1$gene_id$2/ or die "Error, no replacement for identifier $counter, $gene_id";

    }
    
    return($gene_tree_text);
}

