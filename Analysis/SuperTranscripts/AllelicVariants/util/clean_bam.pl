#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\t\tusage: $0 target.bam  errors.txt  output.bam\n\n";

my $target_bam = $ARGV[0] or die $usage;
my $errors_file = $ARGV[1] or die $usage;
my $output_bam = $ARGV[2] or die $usage;



main: {
    
    my %error_reads = &parse_error_reads($errors_file);

    if (%error_reads) {
        
        my $in_cmd = "samtools view -h $target_bam | ";
        
        my $out_cmd = " | samtools view -Sb -o $output_bam";
        
        open(my $in_fh, $in_cmd) or die "Error, cannot open cmd: $in_cmd";
        open(my $out_fh, $out_cmd) or die "Error, cannot open cmd: $out_cmd";
        
        while (my $line = <$in_fh>) {
            if ($line =~ /^\@/) {
                # header line
                print $out_fh $line;
            }
            else {
                my @x = split(/\t/, $line);
                unless ($error_reads{$x[0]}) {
                    print $out_fh $line;
                }
            }
        }
        close $in_fh;
        close $out_fh;
        
        
        &process_cmd("samtools index $output_bam");

        my $output_index = $output_bam;
        $output_index =~ s/\.bam$/\.bai/;

        rename("$output_bam.bai",  $output_index);
        
    }
    else {
        
        &process_cmd("ln -s $target_bam $output_bam");
        
        ###################
        # include the index
        
        $target_bam =~ s/\.bam$/\.bai/;
        $output_bam =~ s/\.bam$/\.bai/;
        
        &process_cmd("ln -s $target_bam $output_bam");

        
        
    }
    
    
    
    exit(0);
}

####
sub process_cmd {
    my ($cmd) = @_;
    
    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, CMD: $cmd died with ret $ret";
    }
    
    return;
}


####
sub parse_error_reads {
    my ($file) = @_;
    
    # examples
    
    #ERROR: Record 6714568, Read name NS500412:143:HJWV5AFXX:1:21307:19673:9975, Read CIGAR M operator maps off end of reference
    #ERROR: Record 6714569, Read name NS500412:143:HJWV5AFXX:3:21405:18461:10525, Read CIGAR M operator maps off end of reference
    #ERROR: Record 6714570, Read name NS500412:143:HJWV5AFXX:1:21307:21758:18659, Alignment start (262144) must be <= reference sequence length (241) on reference TRINITY_DN11029_c1_g1
    #ERROR: Record 6714570, Read name NS500412:143:HJWV5AFXX:1:21307:21758:18659, Read CIGAR M operator maps off end of reference
    #ERROR: Record 6714571, Read name NS500412:143:HJWV5AFXX:2:11202:8325:4497, Alignment start (262144) must be <= reference sequence length (241) on reference TRINITY_DN11029_c1_g1
    
    
    my %error_reads;
    
    open(my $fh, $file) or die "Error, cannot open file: $file";
    while (<$fh>) {
        if (/ERROR: Record (\d+), Read name (\S+), /) {
            my $read_name = $2;
            $error_reads{$read_name} = 1;
        }
    }
    close $fh;
    
    return(%error_reads);
}





    
