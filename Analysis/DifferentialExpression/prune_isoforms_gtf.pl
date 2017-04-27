#!/usr/bin/env perl

use strict;
use warnings;

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use FindBin;
use Data::Dumper;

my $help_flag;


my $usage = <<__EOUSAGE__;

####################################################################################
#
# required:
#
#  --gtf <string>                        gtf transcript structure file
#
#  --pruned_transcripts_fasta <string>   fasta file containing transcript sequences.      
#
####################################################################################

__EOUSAGE__

    ;



my $gtf_file;
my $transcripts_fasta_file;

&GetOptions ( 'h' => \$help_flag,

              'gtf=s' => \$gtf_file,
              'pruned_transcripts_fasta=s' => \$transcripts_fasta_file,
                            
    );


unless ($gtf_file && $transcripts_fasta_file) {
    die $usage;
}



main: {


    my %trans_ids_want;
    {
        open (my $fh, $transcripts_fasta_file) or die "Error, cannot open file $transcripts_fasta_file";
        while (<$fh>) {
            if (/^>(\S+)/) {
                my $trans_id = $1;
                $trans_ids_want{$trans_id} = 1;
            }
        }
        close $fh;
    }

    my %remain_to_find_in_gtf = %trans_ids_want;

    open(my $fh, $gtf_file) or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; }

        my $line = $_;

        if (/transcript_id \"([^\"]+)\"/) {
            my $trans_id = $1;
            
            if ($trans_ids_want{$trans_id}) {
                print $line;
                
                delete $remain_to_find_in_gtf{$trans_id};
                
            }
        }
        
    }
    close $fh;
    
    
    if (%remain_to_find_in_gtf) {

        die "Error, didn't encounter entries for the following transcripts in the gtf file: " . Dumper(\%remain_to_find_in_gtf);
    
    }
    
    exit(0);
    
}

