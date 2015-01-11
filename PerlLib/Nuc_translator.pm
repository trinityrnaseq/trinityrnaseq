#!/usr/bin/env perl

package main;
our $SEE;

package Nuc_translator;

use strict;
require Exporter;
use Carp;


our @ISA = qw (Exporter);
our @EXPORT = qw (translate_sequence get_protein reverse_complement);

use vars qw ($currentCode %codon_table $init_codon_table_subref);


=head1 NAME

package Nuc_translator.pm


=head1 SYNOPSIS

Nuc_translator::use_specified_genetic_code ("universal"); 

my $nuc_sequence = "atgaaagggccctga";

my $translation_frame = 1;

my $protein = &translate_sequence($nuc_sequence, $translation_frame);


=head1 DESCRIPTION

Methods are provided to translate nucleotide sequences into protein sequences using a specified genetic code.

Available genetic codes include universal, Euplotes, Tetrahymena, Candida, Acetabularia

For info on these codes, visit:

http://golgi.harvard.edu/biolinks/gencode.html

Methods exported by this package include:

translate_sequence() 

get_protein() 

reverse_complement()

To change the translation code, the following fully qualified method must be used:

Nuc_translator::use_specified_genetic_code()


=head1 Methods


=cut



## See http://golgi.harvard.edu/biolinks/gencode.html
my %SUPPORTED_GENETIC_CODES = ( universal => 1,
                                Euplotes => 1,
                                Tetrahymena => 1,
                                Candida => 1,
                                Acetabularia => 1,
                                'Mitochondrial-Canonical' => 1,
                                'Mitochondrial-Vertebrates' => 1,
                                'Mitochondrial-Arthropods' => 1,
                                'Mitochondrial-Echinoderms' => 1,
                                'Mitochondrial-Molluscs' => 1,
                                'Mitochondrial-Ascidians' => 1,
                                'Mitochondrial-Nematodes' => 1,
                                'Mitochondrial-Platyhelminths' => 1,
                                'Mitochondrial-Yeasts' => 1,
                                'Mitochondrial-Euascomycetes' => 1,
                                'Mitochondrial-Protozoans' => 1,
                                );



=over 4

=item translate_sequence()

B<Description:> translates a nucleotide sequence given a specific frame 1-6.

B<Parameters:> $nuc_sequence, $frame

B<Returns:> $protein_sequence

=back

=cut



sub translate_sequence {
  my ($sequence, $frame) = @_;
    
    $sequence = uc ($sequence);
	$sequence =~ tr/U/T/;
    my $seq_length = length ($sequence);
    unless ($frame >= 1 and $frame <= 6) { 
		confess "Frame $frame is not allowed. Only between 1 and 6";
	}
	
	if ($frame > 3) {
		# on reverse strand. Revcomp the sequence and reset the frame
		$sequence = &reverse_complement($sequence);
		if ($frame == 4) {
			$frame = 1;
		}
		elsif ($frame == 5) {
			$frame = 2;
		}
		elsif ($frame == 6) {
			$frame = 3;
		}
	}
	
    $sequence =~ tr/T/U/;
    my $start_point = $frame - 1;
    my $protein_sequence;
    for (my $i = $start_point; $i < $seq_length; $i+=3) {
        my $codon = substr($sequence, $i, 3);
        my $amino_acid;
        if (exists($codon_table{$codon})) {
            $amino_acid = $codon_table{$codon};
        } else {
            if (length($codon) == 3) {
                $amino_acid = 'X';
            } else {
                $amino_acid = "";
            }
        }
        $protein_sequence .= $amino_acid;
    }
    return($protein_sequence);
}



=over 4

=item get_protein()

B<Description:> translates nucleotide sequence into a protein sequence.  All 3 forward translation frames are tried
and the first reading frame found to translate without stop codons is returned.  If all 3 frames provide stop codons, the protein with the least number of stops is returned.

B<Parameters:> $nucleotide_sequence

B<Returns:> $protein_sequence

=back

=cut



sub get_protein {
    my ($sequence) = @_;
    
    ## Assume frame 1 unless multiple stops appear.
    my $least_stops = undef();
    my $least_stop_prot_seq = "";
    foreach my $forward_frame (1, 2, 3) {
        my $protein = &translate_sequence($sequence, $forward_frame);
        my $num_stops = &count_stops_in_prot_seq($protein);
        if ($num_stops == 0) {
            return ($protein);
        } else {
            if (!defined($least_stops)) {
                #initialize data
                $least_stops = $num_stops;
                $least_stop_prot_seq = $protein;
            } elsif ($num_stops < $least_stops) {
                $least_stops = $num_stops;
                $least_stop_prot_seq = $protein;
            } else {
                #keeping original $num_stops and $least_stop_prot_seq
            }
        }
    }
    return ($least_stop_prot_seq);
}


=over 4

=item reverse_complement()

B<Description:> reverse complements a nucleotide sequence

B<Parameters:> $nucleotide_sequence

B<Returns:> $nucleotide_sequence_rev_comped

=back

=cut



sub reverse_complement {
    my($s) = @_;
    my ($rc);
    $rc = reverse ($s);
    $rc =~tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
    return($rc);
}


####
sub count_stops_in_prot_seq {
    my ($prot_seq) = @_;
    chop $prot_seq; #remove trailing stop.
    my $stop_num = 0;
    while ($prot_seq =~ /\*/g) {
        $stop_num++;
    } 
    return ($stop_num);
}



####
sub use_specified_genetic_code {
    
    my ($special_code) = @_;
    print STDERR "using special genetic code $special_code\n" if $SEE;
    unless ($SUPPORTED_GENETIC_CODES{$special_code}) {
        die "Sorry, $special_code is not currently supported or recognized.\n";
    }
    &$init_codon_table_subref(); ## Restore default universal code.  Others are variations on this.
    $currentCode = $special_code;
    
    if ($special_code eq "Euplotes") {
        $codon_table{UGA} = "C";
    } 
    
    elsif ($special_code eq "Tetrahymena" || $special_code eq "Acetabularia") {
        $codon_table{UAA} = "Q";
        $codon_table{UAG} = "Q";
    }
    
    elsif ($special_code eq "Candida") {
        $codon_table{CUG} = "S";
    }
    
    elsif ($special_code =~ /Mitochondrial/) {
        &_set_mitochondrial_code($special_code);
    }
 
    else {
        ## shouldn't ever get here anyway.
        confess "Error, code $special_code is not recognized.\n";
    }
    
    
}


####
sub _set_mitochondrial_code {
    my $code = shift;
    ## set canonical by default:
    $codon_table{AUA} = "I";
    $codon_table{AAA} = "K";
    $codon_table{AGA} = $codon_table{AGG} = "R";
    $codon_table{CAU} = $codon_table{CAG} = $codon_table{CAC} = $codon_table{CAA} = "L";
    

    if ($code eq "Mitochondrial-Vertebrates") {
        $codon_table{UGA} = "W";
        $codon_table{AUA} = "M";
        $codon_table{AGA} = "*";
        $codon_table{AGG} = "*";
    }
    elsif ($code eq "Mitochondrial-Arthropods") {
        $codon_table{UGA} = "W";
        $codon_table{AUA} = "M";
        $codon_table{AGA} = "S";
    }
    
    else {
        confess "Sorry, $code hasn't been fully implemented yet.";
    }
   

    ## need to finish


    return;
}



####
sub get_stop_codons {
    my @stop_codons;
    foreach my $codon (keys %codon_table) {
        if ($codon_table{$codon} eq '*') {
            push (@stop_codons, $codon);
        }
    }
    foreach my $codon (@stop_codons) {
        $codon =~ tr/U/T/;
    }
    return (@stop_codons);
}


BEGIN {
  $init_codon_table_subref = sub {
	print STDERR "initing codon table.\n" if $SEE;
    ## Set to Universal Genetic Code
    $currentCode = "universal";
    
    %codon_table = (    UUU => 'F',
                        UUC => 'F',
                        UUA => 'L',
                        UUG => 'L',
                        
                        CUU => 'L',
                        CUC => 'L',
                        CUA => 'L',
                        CUG => 'L',
                        
                        AUU => 'I',
                        AUC => 'I',
                        AUA => 'I',
                        AUG => 'M',
                        
                        GUU => 'V',
                        GUC => 'V',
                        GUA => 'V',
                        GUG => 'V',
                        
                        UCU => 'S',
                        UCC => 'S',
                        UCA => 'S',
                        UCG => 'S',
                        
                        CCU => 'P',
                        CCC => 'P',
                        CCA => 'P',
                        CCG => 'P',
                        
                        ACU => 'T',
                        ACC => 'T',
                        ACA => 'T',
                        ACG => 'T',
                        
                        GCU => 'A',
                        GCC => 'A',
                        GCA => 'A',
                        GCG => 'A',
                        
                        UAU => 'Y',
                        UAC => 'Y',
                        UAA => '*',
                        UAG => '*',
                        
                        CAU => 'H',
                        CAC => 'H',
                        CAA => 'Q',
                        CAG => 'Q',
                        
                        AAU => 'N',
                        AAC => 'N',
                        AAA => 'K',
                        AAG => 'K',
                        
                        GAU => 'D',
                        GAC => 'D',
                        GAA => 'E',
                        GAG => 'E',
                        
                        UGU => 'C',
                        UGC => 'C',
                        UGA => '*',
                        UGG => 'W',
                        
                        CGU => 'R',
                        CGC => 'R',
                        CGA => 'R',
                        CGG => 'R',
                        
                        AGU => 'S',
                        AGC => 'S',
                        AGA => 'R',
                        AGG => 'R',
                        
                        GGU => 'G',
                        GGC => 'G',
                        GGA => 'G',
                        GGG => 'G'    
                        
                        );
  };

  &$init_codon_table_subref();
}


1; #end of module




