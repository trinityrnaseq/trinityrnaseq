#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;
use GFF3_utils;
use SAM_reader;
use SAM_entry;
use Fasta_reader;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

###################################################################
#
# --encoded_fa <string>             feature-encoded fasta file
#
# --coord_sorted_sam <string>       sam alignment file
#
#  Optional:
#
# --SS_lib_type <string>            [FR,RF] (paired), or [R,F] (single)
#
##################################################################

_EOUSAGE_

    ;

my $encoded_fa;
my $sam_file;
my $SS_lib_type;
my $help;

&GetOptions ( 'h' => \$help,
              
              'encoded_fa=s' => \$encoded_fa,
              'coord_sorted_sam=s' => \$sam_file,
              
              'SS_lib_type=s' => \$SS_lib_type,
              );


if ($help) {
    die $usage;
}

if (! ($encoded_fa && $sam_file) ) {
    die $usage;
}



my %SS_encoding = ( 0 => 'intergenic',
                    1 => 'intron',
                    2 => 'exon+',
                    3 => 'exon-',
                    4 => 'rRNA+',
                    5 => 'rRNA-',
                    );

my %reg_encoding = ( 0 => 'intergenic',
                     1 => 'intron',
                     2 => 'exon',
                     3 => 'exon',
                     4 => 'rRNA',
                     5 => 'rRNA',
                     );



main: {


    print STDERR "-parsing genome encoding\n";
    my $fasta_reader = new Fasta_reader($encoded_fa);
    my %genome_encoding = $fasta_reader->retrieve_all_seqs_hash();
    
    my %feature_coverage_counter;


    print STDERR "-parsing sam file\n";
    
    my $curr_scaffold = "";
    my $feat_encoding = "";
    #my @feats;
    
    my $sam_reader = new SAM_reader($sam_file);
    
    my $read_counter = 0;
    while (my $sam_entry = $sam_reader->get_next()) {
        
        $read_counter++;
        print STDERR "\r[$read_counter] reads processed.    " if ($read_counter % 1000 == 0);

        if ($sam_entry->is_query_unmapped()) {
            next;
        }
        

        my $scaffold = $sam_entry->get_scaffold_name();
        if ($scaffold ne $curr_scaffold) {
            $curr_scaffold = $scaffold;
            $feat_encoding = $genome_encoding{$curr_scaffold} || ""; # no features on a contig are considered intergenic
            #@feats = split(//, $feat_encoding);
            #print "encoding: $feat_encoding, length: " . length($feat_encoding) . "\n";
        }
        
        my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();
       
        my $strand = ($SS_lib_type) ? $sam_entry->get_query_transcribed_strand($SS_lib_type) : $sam_entry->get_query_strand();
        
        my $genome_coordset = shift @$genome_coords_aref;
        my ($lend, $rend) = @$genome_coordset;
        my $midpt = int( ($rend+$lend)/2);
        my $encoding = 0;
        if ($midpt -1 < length($feat_encoding)) {
            $encoding = substr($feat_encoding, $midpt-1, 1);
        }
        
        my $feature_type = ($SS_lib_type) ? $SS_encoding{$encoding} : $reg_encoding{$encoding};
        
        if ($SS_lib_type) {
            if ($feature_type =~ /([\+\-])$/) {
                my $feature_orient = $1;
                if ($feature_orient eq $strand && $strand eq '-') {
                    #  -,-
                    # consider it a forward feature
                    $feature_type =~ s/\-$/\+/;
                }
                elsif ($feature_orient ne $strand && $feature_orient eq '+') {
                    # Feature+,read-  : make feature -
                    $feature_type =~ s/\+$/-/;
                }
                # Feature-, read+  : counting feature as anti, so already OK
                # Feature+, read+  : counting feature as plus, so already OK
            }
            $feature_coverage_counter{$feature_type}++;
        }
        else {
            ## not strand-specific
            $feature_coverage_counter{$feature_type}++;
        }
    }
    
    my $total_coverage = 0;
    foreach my $count (values %feature_coverage_counter) {
        $total_coverage += $count;
    }


    print "\n\n";
    foreach my $feature_type (sort keys %feature_coverage_counter) {
        my $coverage = $feature_coverage_counter{$feature_type};
        my $percent = sprintf("%.2f", $coverage/$total_coverage*100);
        
        print join("\t", $feature_type, $coverage, "$percent\%") . "\n";
    }
    print "\n\n";
    
    exit(0);
}

