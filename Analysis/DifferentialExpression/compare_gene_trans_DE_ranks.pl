#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

my $help_flag;

my $usage = <<__EOUSAGE__;

#############################################################################################
#
#  --genes|G <string>        gene DE_results file
#  --trans|T <string>        trans DE_results file
#
#  --max_fdr|F <float>       maximum FDR (ignore those that dont meet this requirement) (ie. 0.001)
#
#  --max_count_select|M <int>    number of top ranked trans and genes to capture and report. (ie. 1000)
#
#############################################################################################


__EOUSAGE__

    ;




my $MAX_COUNT;
my $MAX_FDR;
my $genes_file;
my $trans_file;

&GetOptions ( 'h' => \$help_flag,
              'genes|G=s' => \$genes_file,
              'trans|T=s' => \$trans_file,
              'max_fdr|F=f' => \$MAX_FDR,
              'max_count_select|M=i' => \$MAX_COUNT,
    );



if ($help_flag) {
    die $usage;
}

unless ($genes_file && $trans_file && $MAX_FDR && $MAX_COUNT) {
    die $usage;
}



my %data;




main: {
    
    &add_info($genes_file, 'gene');
    &add_info($trans_file, 'trans');
    
    print "\tgene_fdr\ttrans_fdr\tgene_logfc\ttrans_logfc\n";
    
    my %gene_to_pt;

    my @pts;
    foreach my $acc (keys %data) {
        my $gene_fdr = $data{$acc}->{gene}->{fdr};
        my $trans_fdr = $data{$acc}->{trans}->{fdr};

        my $gene_logfc = $data{$acc}->{gene}->{logfc} || "NA";
        my $trans_logfc = $data{$acc}->{trans}->{logfc} || "NA";
        

        if (  (defined($gene_fdr) && $gene_fdr <= $MAX_FDR)
              ||
              (defined($trans_fdr) && $trans_fdr <= $MAX_FDR) ) {

            $gene_fdr = (defined $gene_fdr) ? &make_logscale($gene_fdr) : -1;
            $trans_fdr = (defined $trans_fdr) ? &make_logscale($trans_fdr) : -1;
            
            my $sum = 0;
            $sum += $gene_fdr;
            $sum += $trans_fdr;
            
            my $struct = { 

                acc => $acc,
                
                gene_fdr => $gene_fdr,
                trans_fdr => $trans_fdr,
                sum_fdr => $sum,

                gene_logfc => $gene_logfc,
                trans_logfc => $trans_logfc,
            };
            

            $gene_to_pt{$acc} = $struct;

            push (@pts, $struct); 
        }
    }

    
    my $max_pts = scalar(keys %data);
    if ($max_pts < $MAX_COUNT) {
        $MAX_COUNT = $max_pts;
    }
    
    my %top_entries;
    
    # sort by gene
    my $gene_count_by_gene = 0;
    @pts = reverse sort {$a->{gene_fdr}<=>$b->{gene_fdr}} @pts;
    

    use Data::Dumper;
    #print STDERR Dumper(\@pts);
    #die;

    for my $pt (@pts) {
        my $acc = $pt->{acc};
        
        my $gene_fdr = $pt->{gene_fdr};
        
        if ($gene_fdr > 0) {
            $gene_count_by_gene++;
        
            if ($gene_count_by_gene < $MAX_COUNT) {
                $top_entries{$acc}++;
                #print STDERR "$acc => gene\n";
            }
        }
        
    }
    
    # sort by trans
    my $gene_count_by_trans = 0;
    @pts = reverse sort {$a->{trans_fdr}<=>$b->{trans_fdr}} @pts;
    for my $pt (@pts) {
        my $acc = $pt->{acc};
        if ($pt->{trans_fdr} > 0) {
           
            $gene_count_by_trans++;
            
            if ($gene_count_by_trans < $MAX_COUNT) {
                 $top_entries{$acc}++;
                 #print STDERR "$acc => trans\n";
            }
        }
    }
    
    

    my @top_structs;
    my $count_both = 0;
    my $count_either = 0;
    foreach my $acc (keys %top_entries) {
        my $struct = $gene_to_pt{$acc};
        push (@top_structs, $struct);
        
        if ($top_entries{$acc} > 1) {
            $count_both++;
        }
        $count_either++;
    }
    
    my $percent_top_rank_overlap = sprintf("%.2f", $count_both/$count_either * 100);
    print STDERR "counts of genes by genes <= FDR=$MAX_FDR = $gene_count_by_gene\n";
    print STDERR "counts of genes by trans <= FDR=$MAX_FDR = $gene_count_by_trans\n";
    print STDERR "% top rank overlap of top $MAX_COUNT entries = $count_both/$count_either * 100 = $percent_top_rank_overlap %\n";
    
    @top_structs = reverse sort {$a->{sum_fdr}<=>$b->{sum_fdr}} @top_structs;
        
    foreach my $struct (@top_structs) {
        
        my ($acc, $gene_fdr, $trans_fdr, $gene_logfc, $trans_logfc) = ($struct->{acc},
                                                                       $struct->{gene_fdr},
                                                                       $struct->{trans_fdr},
                                                                       $struct->{gene_logfc},
                                                                       $struct->{trans_logfc});
        

        if ($gene_fdr < 0) { $gene_fdr = "NA"; }
        if ($trans_fdr < 0) { $trans_fdr = "NA"; }
        

        print join("\t", $acc, $gene_fdr, $trans_fdr, $gene_logfc, $trans_logfc) . "\n";
    }
    

    exit(0);
    
   

}

####
sub make_logscale {
    my ($fdr) = @_;
    
    $fdr = -1 * (log($fdr)/log(10));
    
    return($fdr);
}

####
sub add_info {
    my ($file, $type) = @_;

    open (my $fh, $file) or die $!;
    my $header = <$fh>;
    my $count = 0;

    while (<$fh>) {
        my @x = split(/\t/);
        my $acc = $x[0];
        my $logfc = $x[1];
        my $fdr = $x[4];
        
        if ($fdr > $MAX_FDR) { last; }

        if ($acc =~ /(comp\d+_c\d+)/) {
            $acc = $1;
        }
        else {
            die "Error, cannot parse component name from $acc";
        }
        
        if (! exists ($data{$acc}->{$type})) {
            $data{$acc}->{$type}->{fdr} = $fdr;
            $data{$acc}->{$type}->{logfc} = $logfc;
            
        }
        

        
        
        

    }
    
    return;
}
