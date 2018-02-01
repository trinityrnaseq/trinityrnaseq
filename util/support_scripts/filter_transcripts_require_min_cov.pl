#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use DelimParser;
use Carp;
use Data::Dumper;

my $usage = "\n\n\tusage: $0 Trinity.fasta reads.fa salmon.quant.sf min_cov\n\n";

my $trin_fa = $ARGV[0] or die $usage;
my $reads_fa = $ARGV[1] or die $usage;
my $salmon_quant_file = $ARGV[2] or die $usage;
my $min_cov = $ARGV[3] or die $usage;

main: {

    my $read_length = &estimate_read_length($reads_fa);
    
    my %want;
    {
        open(my $fh, "$salmon_quant_file") or die "Error, cannot open file: $salmon_quant_file";
        my $tab_reader = new DelimParser::Reader($fh, "\t");
        my %col_headers = map { + $_ => 1 } $tab_reader->get_column_headers();
        # Name    Length  EffectiveLength TPM     NumReads
        unless (exists $col_headers{'Name'} && exists $col_headers{'NumReads'} && exists $col_headers{'Length'}) {
            confess "ERROR, not recognizing salmon header.  Need Name, NumReads, and Length columns";
        }
        while(my $row = $tab_reader->get_row()) {
            my $transcript = $row->{Name};
            my $length = $row->{Length};
            my $count = $row->{NumReads};

            
            unless (defined $transcript) {
                confess "Error, no transcript name provided at " . Dumper($row);
            }
            unless (defined $count) {
                confess "Error, no NumReads value provided at " . Dumper($row);
            }
            unless (defined $length) {
                confess "Error, no length value provided at " . Dumper($row);
            }

            my $eff_cov = $count * $read_length / $length;
            
            if ($eff_cov >= $min_cov) {
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


####
sub estimate_read_length {
    my ($fa_file) = @_;

    my $fasta_reader = new Fasta_reader($fa_file);

    my $max_records = 100;

    my $record_counter = 0;
    my $sum_lens = 0;
    while (my $seq_obj = $fasta_reader->next()) {
        my $sequence = $seq_obj->get_sequence();
        $sum_lens += length($sequence);

        $record_counter++;
        if ($record_counter >= $max_records) {
            last;
        }
    }

    my $avg_len = $sum_lens / $record_counter;
    
    return($avg_len);
}

        
        
        
        
            
