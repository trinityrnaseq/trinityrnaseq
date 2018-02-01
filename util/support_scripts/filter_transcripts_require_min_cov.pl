#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use DelimParser;
use Carp;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 Trinity.fasta salmon.quant.sf min_reads\n\n";

my $trin_fa = $ARGV[0] or die $usage;
my $salmon_quant_file = $ARGV[1] or die $usage;
my $min_reads = $ARGV[2] or die $usage;

main: {

    my %want;
    {
        open(my $fh, "$salmon_quant_file") or die "Error, cannot open file: $salmon_quant_file";
        my $tab_reader = new DelimParser::Reader($fh, "\t");
        my %col_headers = map { + $_ => 1 } $tab_reader->get_column_headers();
        # Name    Length  EffectiveLength TPM     NumReads
        unless (exists $col_headers{'Name'} && exists $col_headers{'TPM'}) {
            confess "ERROR, not recognizing salmon header.  Need TPM and Name columns";
        }
        while(my $row = $tab_reader->get_row()) {
            my $transcript = $row->{Name};
            my $count = $row->{NumReads};

            unless (defined $transcript) {
                confess "Error, no transcript name provided at " . Dumper($row);
            }
            unless (defined $count) {
                confess "Error, no NumReads value provided at " . Dumper($row);
            }

            if ($count >= $min_reads) {
                $want{$transcript} = 1;
            }
        }
    }

    my $fasta_reader = new Fasta_reader($trin_fa);
    while (my $seq_obj = $fasta_reader->next()) {
        my $acc = $seq_obj->get_accession();
        if ($want{$acc}) {
            my $fasta_entry = $seq_obj->get_FASTA_format();

            print $fasta_entry;

            delete $want{$acc};
        }
        
        
    }

    if (%want) {
        confess "Error, missed retrieving entries from fasta file: " . Dumper(%want);
    }
    
    exit(0);
}


        
        
        
        
            
