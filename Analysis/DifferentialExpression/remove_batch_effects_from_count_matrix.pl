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


my $usage = <<__EOUSAGE__;


#################################################################################################
#
#  Required:
#
#  --matrix|m <string>               matrix of raw read counts (not normalized!)
#
#  --batches_file|b <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        batch_1    cond_A_rep1
#                                        batch_1    cond_B_rep1
#
#                                        batch_2    cond_A_rep2
#                                        batch_2    cond_B_rep2
#
#
################################################################################################



__EOUSAGE__


    ;


my $matrix_file;
my $batches_file;
my $help_flag;

&GetOptions ( 'h' => \$help_flag,
              'matrix|m=s' => \$matrix_file,              
              'batches_file|b=s' => \$batches_file,
              
    );



if ($help_flag) {
    die $usage;
}


unless ($matrix_file 
        && $batches_file
    ) { 
    
    die $usage;
    
}

main: {
            
    &run_edgeR_remove_batch_effect($matrix_file, $batches_file);
    
    
    exit(0);
}


####
sub run_edgeR_remove_batch_effect {
    my ($matrix_file, $batches_file) = @_;

    use File::Basename;
    my $base = basename($batches_file);
    my $Rscript_name = "$base.batch_eff_removal.Rscript";
    my $out_matrix_name = basename($matrix_file) . ".batch_eff_removal.matrix";

    ## write R-script to run edgeR
    open (my $ofh, ">$Rscript_name") or die "Error, cannot write to $Rscript_name";
    
    print $ofh "library(edgeR)\n";

    print $ofh "\n";
    
    print $ofh "rnaseqMatrix = read.table(\"$matrix_file\", header=T, row.names=1, com='', check.names=F)\n";
    print $ofh "batches = read.table(\"$batches_file\", header=F, row.names=2, check.names=F)\n";
    print $ofh "batch_factors = as.factor(batches[colnames(rnaseqMatrix),])\n";
    
    print $ofh "exp_study = DGEList(counts=rnaseqMatrix)\n";
    print $ofh "exp_study = calcNormFactors(exp_study)\n";
    print $ofh "logCPM <- cpm(exp_study,log=TRUE)\n";
    print $ofh "logCPM <- removeBatchEffect(logCPM,batch=batch_factors)\n";
    
    print $ofh "millions = colSums(rnaseqMatrix)/1e6\n";
    
    print $ofh "myCPM = 2^logCPM\n";
    print $ofh "myRecounts = apply(myCPM, 1, function (x) x*millions)\n";
    print $ofh "write.table(t(myRecounts), file=\'$out_matrix_name\', quote=F, sep='\t')\n";
        
    
    close $ofh;

    ## Run R-script
    my $cmd = "R --vanilla -q < $Rscript_name";


    eval {
        &process_cmd($cmd);
    };
    if ($@) {
        print STDERR "$@\n\n";
        print STDERR "\n\nWARNING: This EdgeR comparison failed...\n\n";
        ## if this is due to data paucity, such as in small batch data sets, then ignore for now.
    }
    

    return;
}



####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";
    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret ($ret) ";
    }

    return;
}
