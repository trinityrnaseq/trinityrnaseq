#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);

my $help_flag;


my $usage = <<__EOUSAGE;

#################################################
#
# usage:  $0 [opts] < trinity_fasta_files_listing
#
#################################################
#
# Required
#
#  --token_prefix <string>       eg.  TRINITY_DN or TRINITY_GG
#
#  --output_prefix <string>      eg.  Trinity  or Trinity_GG
#
# Optional:
#
#  --include_supertranscripts    flag, captures the supertranscripts fasta and gtf
#
##################################################  




__EOUSAGE

    ;
    

my $token_prefix;
my $output_prefix;
my $include_supertranscripts_flag;

&GetOptions ( 'h' => \$help_flag,
              'token_prefix=s' => \$token_prefix,
              'output_prefix=s' => \$output_prefix,
              'include_supertranscripts' => \$include_supertranscripts_flag,
    );


if ($help_flag) {
    die $usage;
}

unless ($token_prefix && $output_prefix) {
    die $usage;
}

my $ofh_trinity_fasta;
open($ofh_trinity_fasta, ">$output_prefix.fasta") or die "Error, cannot write to $output_prefix.fasta";

my $ofh_supertrans_fasta;
my $ofh_supertrans_gtf;
if ($include_supertranscripts_flag) {
    open($ofh_supertrans_fasta, ">$output_prefix.SuperTrans.fasta") or die $!;
    open($ofh_supertrans_gtf, ">$output_prefix.SuperTrans.gtf") or die $!;
}


main: {

    while (<STDIN>) {
        my $filename = $_;
        chomp $filename;
        unless (-e $filename) {
            print STDERR "ERROR, filename: $filename is indicated to not exist.\n";
            next;
        }
        if (-s $filename) {
            
            # ie. read_partitions/Fb_0/CBin_161/c16167.trinity.reads.fa
            
            $filename =~ m|c(\d+)\.trinity\.reads| or die "Error, cannot parse Trinity component value from filename: $filename";
            my $component = $1;

            &process_files($filename, $component);
            
        }
    }

    exit(0);
}
            

####
sub process_files {
    my ($filename, $component) = @_;

    {
        # write fasta
        open (my $fh, $filename) or die "Error, cannot open file $filename";
        while (<$fh>) {
            if (/>/) {
                s/>/>${token_prefix}${component}\_/;
            }
            print $ofh_trinity_fasta $_;
        }
        close $fh;
    }

    if ($include_supertranscripts_flag) {
        # write supertranscripts
        
        my $core_filename = $filename;
        $core_filename =~ s/\.fasta$//;

        { # supertranscripts fasta:
            my $supertranscripts_fasta = "$core_filename.SuperTrans.fasta";
            if (! -s $supertranscripts_fasta) {
                die "Error, missing file: $supertranscripts_fasta";
            }
            open(my $fh, $supertranscripts_fasta) or die "Error, cannot open file: $supertranscripts_fasta";
            while (<$fh>) {
                if (/>/) {
                    s/>/>${token_prefix}${component}\_/;
                }
                print $ofh_supertrans_fasta $_;
            }
            close $fh;
        }

        {
            # supertranscripts gtf
            my $supertranscripts_gtf = "$core_filename.SuperTrans.gtf";
            if (! -s $supertranscripts_gtf) {
                die "Error, missing file: $supertranscripts_gtf";
            }
            open(my $fh, $supertranscripts_gtf) or die "Error, cannot open file $supertranscripts_gtf";
            while (<$fh>) {
                my $line = $_;
                if (/\w/) {

                    my @x = split(/\t/, $line);
                    my $id = $x[0];
                    $line =~ s/$id/${token_prefix}${component}_${id}/g;

                }
                print $ofh_supertrans_gtf $line;
            }
            close $fh;
        }
    }
    
    return;
}

