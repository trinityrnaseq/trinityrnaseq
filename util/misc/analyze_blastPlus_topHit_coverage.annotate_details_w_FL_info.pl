#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $help_flag;

my $usage = <<__EOUSAGE__;

#################################################################################################
#
# Required:
#
#  --details <string>              compreh_init_build.details filename
#  --blast <string>                blast.outfmt6.w_pct_hit_length filename
#
# Optional:
#
#  --min_pct_hit_len <float>     min length of alignment coverage of hit (default: 80)
#  --pasa_validations <string>   incorporate notes from pasa's alignment validations file (alignment.validations.out)
#
#################################################################################################


__EOUSAGE__

    ;


my $compreh_details_file;
my $blast_file;
my $min_pct_hit_len = 80;
my $pasa_validations_file;

&GetOptions ( 'h' => \$help_flag,
              'details=s' => \$compreh_details_file,
              'blast=s' => \$blast_file,
              'min_pct_hit_len=f' => \$min_pct_hit_len,
              'pasa_validations=s' => \$pasa_validations_file,
              
              );


if ($help_flag) {
    die $usage;
}

unless ($compreh_details_file && $blast_file) {
    die $usage;
}

my %pasa_validations_info;
if ($pasa_validations_file) {
    
    open (my $fh, $pasa_validations_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $acc = $x[1];
        my $note = $x[14];
        if ($note) {
            $pasa_validations_info{$acc} .= $note . "; ";
        }
    }
    close $fh;

}

my %FL_mappings;

{
    open (my $fh, $blast_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        
        chomp;
        my @x = split(/\t/);
        my $pct_hit_len = $x[13];
        if ($pct_hit_len < $min_pct_hit_len) {
            next;
        }
        
        my $query = $x[0];
        my $hit = $x[1];
        my $annot = $x[14];
        my $Evalue = $x[10];

        $FL_mappings{$query} = join("\t", $hit, $annot, $Evalue, $pct_hit_len);
        
    }
    close $fh;
}

my $prev_gene = "";
open (my $fh, $compreh_details_file) or die $!;
while (<$fh>) {
    chomp;
    my @x = split(/\t/);
    my $gene = $x[0];
    my $acc = $x[1];
    if (my $mappings = $FL_mappings{$acc}) {
        $mappings =~ s/\t/ /g;
        push (@x, $mappings);
    }
    else {
        push (@x, "");
    }
    if (my $validation_note = $pasa_validations_info{$acc}) {
        push (@x, "validation_note: $validation_note");
    }
    

    if ($gene ne $prev_gene) {
        print "\n";
    }
    $prev_gene = $gene;
    
    print join("\t", @x) . "\n";
}



exit(0);


