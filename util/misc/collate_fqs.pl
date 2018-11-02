#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

use FindBin;
use lib("$FindBin::Bin/../../PerlLib");
use Fastq_reader;
use Data::Dumper;

my $help_flag;

my $usage = <<__EOUSAGE__;

#############################################################################################
#
# --samples_file <string>      samples.txt
#
# --output_prefix <string>     outputs will be named <string>.left.fq and <string>.right.fq
#
#############################################################################################



__EOUSAGE__

    ;


my $samples_file;
my $output_prefix;

&GetOptions ( 'h' => \$help_flag,
              'samples_file=s' => \$samples_file,
              'output_prefix=s' => \$output_prefix,
    );

unless ($samples_file && $output_prefix) {
    die $usage;
}

main: {
    
    my @samples = &parse_samples_file($samples_file);

    # init fq readers    
    foreach my $sample (@samples) {
        my $left_fq = $sample->{left_fq};
        my $left_fq_reader = new Fastq_reader($left_fq);
        $sample->{left_fq_reader} = $left_fq_reader;

        my $right_fq = $sample->{right_fq};
        my $right_fq_reader = new Fastq_reader($right_fq);
        $sample->{right_fq_reader} = $right_fq_reader;
    }

    # open output files:
    my $left_fq_outfile = "$output_prefix.left.fq";
    open(my $left_ofh, ">$left_fq_outfile") or die "Error, cannot write to $left_fq_outfile";

    my $right_fq_outfile = "$output_prefix.right.fq";
    open(my $right_ofh, ">$right_fq_outfile") or die "Error, cannot write to $right_fq_outfile";
    
    my $samples_remaining = scalar(@samples);

    my $counter = 0;
    while($samples_remaining != 0) {
                
        $samples_remaining = 0;
        
        foreach my $sample (@samples) {

            my $sample_name = $sample->{sample_name};
            
            my $left_fq_reader = $sample->{left_fq_reader};
            my $left_fq_entry = $left_fq_reader->next();

            my $right_fq_reader = $sample->{right_fq_reader};
            my $right_fq_entry = $right_fq_reader->next();

            if ($left_fq_entry xor $right_fq_entry) {
                confess "Error, found left_fq entry but not right_fq entry: " . Dumper([$sample, $left_fq_entry, $right_fq_entry]);
            }
            elsif ($left_fq_entry && $right_fq_entry) {
                                
                print $left_ofh $left_fq_entry->get_fastq_record();
                print $right_ofh $right_fq_entry->get_fastq_record();
                
                $samples_remaining = 1;
            }
            
        }

        $counter++;

        if ($counter % 100000 == 0) {
            print STDERR "\r[$counter]        ";
        }
    } # end of while samples remaining

    close $left_ofh;
    close $right_ofh;

    print STDERR "\n\nDone.\n";
    
    exit(0);
    
}




####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;
    
    open(my $fh, $samples_file) or die "Error, cannot open file: $samples_file";
    while(<$fh>) {
        chomp;
        my ($sample_name, $left_fq, $right_fq) = split(/\t/);
        push (@samples, { sample_name => $sample_name,
                          left_fq => $left_fq,
                          right_fq => $right_fq,
              } );
    }

    close $fh;

    return(@samples);
}

