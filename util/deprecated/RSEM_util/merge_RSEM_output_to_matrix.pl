#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);

use File::Basename;

my $usage = <<__EOUSAGE__;


############################################################################################
#
# Required:
#
#  --rsem_files <string>       file containing a list of RSEM output files.
#                                  (note, should be isoform or gene.results files, dont mix)
#
#  --mode <string>             counts|fpkm|tpm
#
# Optional:
#
#  --trans_mode_bundle_gene_id    specify the transcript identifier as 'gene|trans'
#
############################################################################################

__EOUSAGE__

    ;



my $rsem_files_list_file;
my $mode;
my $trans_mode_bundle_gene_id = 0;
my $help_flag;

&GetOptions( 'help|h' => \$help_flag,
             
             'rsem_files=s' => \$rsem_files_list_file,
             'mode=s' => \$mode,
             'trans_mode_bundle_gene_id' => \$trans_mode_bundle_gene_id,
             );

if ($help_flag) {
    die $usage;
}
unless ($rsem_files_list_file && $mode) {

    die $usage;
}

unless ($mode =~ /^(counts|fpkm|tpm)$/i) {
    die "Error, mode $mode not recognized";
}
                                             

=header_format

0       transcript_id
1       gene_id
2       length
3       effective_length
4       expected_count
5       TPM
6       FPKM
7       IsoPct

=cut


main: {


    unless (-f $rsem_files_list_file) { die "$usage\nError: cannot open file $rsem_files_list_file"; }
    my @rsem_files = `cat $rsem_files_list_file`;
    chomp @rsem_files;
    
    my %data;
    
    foreach my $file (@rsem_files) {
        
        print STDERR "-capturing $mode from file: $file\n";
        
        open (my $fh, $file) or die "Error, cannot open file $file";
        my $header = <$fh>; # ignore it
        while (<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];

            if ($trans_mode_bundle_gene_id) {
                my $gene = $x[1];
                $acc = "$gene|$acc";
            }
            
            my $count = $x[4];
            my $tpm = $x[5];
            my $fpkm = $x[6];

            
            $data{$acc}->{$file} = ($mode =~ /counts/i) ? $count : ($mode =~ /fpkm/) ? $fpkm : $tpm;
        }
        close $fh;
    }

    
    print STDERR "\n\n* done parsing files.  Now outputting matrix.\n\n";
    
    my @filenames = @rsem_files;
    foreach my $file (@filenames) {
        $file = basename($file);
        $file =~ s/-/_/g; # R doesn't like '-' in column headers
        $file =~ s/\.(genes|isoforms)\.results//;
    }
    
    
    print join("\t", "", @filenames) . "\n";
    foreach my $acc (keys %data) {
        
        print "$acc";

        foreach my $file (@rsem_files) {

            my $count = $data{$acc}->{$file};
            unless (defined $count) {
                $count = "NA";
            }

            print "\t$count";
            
        }
        
        print "\n";
        
    }
    
    
    exit(0);
}
    
        
