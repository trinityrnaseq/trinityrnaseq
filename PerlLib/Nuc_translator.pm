#!/usr/bin/env perl

package main;
our $SEE;

package Nuc_translator;

use strict;
require Exporter;
use Carp;


our @ISA = qw (Exporter);
our @EXPORT = qw (get_genetic_codes translate_sequence get_protein reverse_complement);

use vars qw ($currentCode %codon_table $init_codon_table_subref $support_Thymine_and_Uracil_subref);


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
(https://web.archive.org/web/20040216020103/http://golgi.harvard.edu/biolinks/gencode.html)

Methods exported by this package include:

translate_sequence() 

get_protein() 

reverse_complement()

To change the translation code, the following fully qualified method must be used:

Nuc_translator::use_specified_genetic_code()


=head1 Methods


=cut



## See http://golgi.harvard.edu/biolinks/gencode.html
my %SUPPORTED_GENETIC_CODES = ( Universal => 1,
                                
                                Tetrahymena => 1,
                                Acetabularia => 1,
                                Ciliate => 1,
                                Dasycladacean => 1,
                                Hexamita => 1,
                                
                                Candida => 1,
                                
                                Euplotid => 1,

                                SR1_Gracilibacteria => 1,

                                Pachysolen_tannophilus => 1,

                                Mesodinium => 1,

                                Peritrich => 1,
                                
                                
                                
                                'Mitochondrial-Vertebrates' => 1,
                                'Mitochondrial-Yeast' => 1,
                                'Mitochondrial-Invertebrates' => 1,
                                "Mitochondrial-Protozoan" => 1,
                                "Mitochondrial-Echinoderm" => 1,
                                "Mitochondrial-Ascidian" => 1,
                                "Mitochondrial-Flatworm" => 1,
                                "Mitochondrial-Chlorophycean" => 1,
                                "Mitochondrial-Trematode" => 1,
                                "Mitochondrial-Scenedesmus_obliquus" => 1,
                                "Mitochondrial-Thraustochytrium" => 1,
                                "Mitochondrial-Pterobranchia" => 1,
                                
                                );





=over 4

=item get_genetic_codes()

B<Description:> provides the list of supported genetic codes

B<Parameters:> 

B<Returns:> list of genetic codes

=back

=cut

####
sub get_genetic_codes {
    return (sort keys %SUPPORTED_GENETIC_CODES);
}



####
sub show_translation_table {
    my ($genetic_code) = @_;
    
    &use_specified_genetic_code($genetic_code);

    my @nucs = qw(T C A G);
    

    my $translation_line = "";
    my $c1_line = "";
    my $c2_line = "";
    my $c3_line = "";
    
    for my $c1 (@nucs) {
        for my $c2 (@nucs) {
            for my $c3 (@nucs) {
                $c1_line .= $c1;
                $c2_line .= $c2;
                $c3_line .= $c3;

                my $codon = $c1 . $c2 . $c3;
                my $translation = $codon_table{$codon};
                $translation_line .= "$translation";
            }
        }
    }

    my $translation_table = join("\n", $translation_line, $c1_line, $c2_line, $c3_line) . "\n";

    return($translation_table);
}


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
    $rc =~tr/ACGTacgtyrkmYRKMUu/TGCAtgcarymkRYMKAa/;
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

    if ($special_code eq "Universal") {
        # already set.
    }
        
    elsif ($special_code =~ /^(Tetrahymena|Acetabularia|Ciliate|Dasycladacean|Hexamita)$/) {
        $codon_table{UAA} = "Q";
        $codon_table{UAG} = "Q";
    }
    
    elsif ($special_code eq "Candida") {
        $codon_table{CUG} = "S";
    }

    elsif ($special_code eq "Euplotid") {
        $codon_table{UGA} = "C"; # not *
    }

    elsif ($special_code eq "SR1_Gracilibacteria") {
        $codon_table{UGA} = "G"; # not *
    }

    elsif ($special_code eq "Pachysolen_tannophilus") {
        $codon_table{CUG} = "A"; # not L
    }
    elsif ($special_code eq "Mesodinium") {
        $codon_table{UAA} = "Y"; # not *
        $codon_table{UAG} = "Y"; # not *
    }

    elsif ($special_code eq "Peritrich") {
        $codon_table{UAA} = "E"; # not *
        $codon_table{UAG} = "E"; # not *
    }
    
    elsif ($special_code =~ /Mitochondrial/) {
        &_set_mitochondrial_code($special_code);
    }
    
    else {
        ## shouldn't ever get here anyway.
        confess "Error, code $special_code is not recognized.\n";
    }
    

    &$support_Thymine_and_Uracil_subref();
    
}


####
sub _set_mitochondrial_code {
    my $code = shift;
        
    # see: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2
    
    if ($code eq "Mitochondrial-Vertebrates") {
        $codon_table{UGA} = "W";
        $codon_table{AUA} = "M";
        $codon_table{AGA} = "*";
        $codon_table{AGG} = "*";
    }
    elsif ($code eq "Mitochondrial-Yeast") {
        $codon_table{AUA} = "M"; # instead of I
        $codon_table{CUU} = $codon_table{CUC} = $codon_table{CUA} = $codon_table{CUG} = "T"; # instead of L
        $codon_table{UGA} = "W"; # instead of *
    }
    elsif ($code eq "Mitochondrial-Invertebrates") {
        $codon_table{AGA} = "S";
        $codon_table{AGG} = "S";
        $codon_table{AUA} = "M";        
        $codon_table{UGA} = "W";
    }
    elsif ($code eq "Mitochondrial-Protozoan") {
        $codon_table{UGA} = "W";
    }
    elsif ($code eq "Mitochondrial-Echinoderm" || $code eq "Mitochondrial-Flatworm") {
        $codon_table{AAA} = "N"; # not K
        $codon_table{AGA} = "S"; # not R
        $codon_table{AGG} = "S"; # not R
        $codon_table{UGA} = "W"; # not *
    }
    elsif ($code eq "Mitochondrial-Ascidian") {
        $codon_table{AGA} = "G"; # not R
        $codon_table{AGG} = "G"; # not R
        $codon_table{AUA} = "M"; # not I
        $codon_table{UGA} = "W"; # not *
    }
    
    elsif ($code eq "Mitochondrial-Chlorophycean") {
        $codon_table{UAG} = "L"; # not *
    }
    elsif ($code eq "Mitochondrial-Trematode") {
        $codon_table{UGA} = "W"; # not *
        $codon_table{AUA} = "M"; # not I
        $codon_table{AGA} = "S"; # not R
        $codon_table{AGG} = "S"; # not R
        $codon_table{AAA} = "N"; # not K
    }
    elsif ($code eq "Mitochondrial-Scenedesmus_obliquus") {
        $codon_table{UCA} = "*"; # not S
        $codon_table{UAG} = "L"; # not *
    }
    elsif ($code eq "Mitochondrial-Thraustochytrium") {
        $codon_table{UUA} = "*";
    }
    elsif ($code eq "Mitochondrial-Pterobranchia") {
        $codon_table{AGA} = "S"; # not R
        $codon_table{AGG} = "K"; # not R
        $codon_table{UGA} = "W"; # not *
    }
    
    
    else {
        confess "Sorry, $code hasn't been fully implemented yet.";
    }
    

    

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
    #foreach my $codon (@stop_codons) {
    #    $codon =~ tr/U/T/;
    #}
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


  $support_Thymine_and_Uracil_subref = sub {

      my @codons = keys %codon_table;
      foreach my $codon (@codons) {
          if ($codon =~ /U/) {
              my $aa = $codon_table{$codon};
              my $T_codon = $codon;
              $T_codon =~ s/U/T/g;
              $codon_table{$T_codon} = $aa;
          }
      }
  };
  
  
  # init codon table, using uracil codons
  &$init_codon_table_subref();

  # update codon table to also support thymine
  &$support_Thymine_and_Uracil_subref();
}




####
sub run_test() {
    my %expected_translation_tables = (
        'Universal' => "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Vertebrates' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        'Mitochondrial-Yeast' => "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Protozoan' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Invertebrates' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        'Ciliate' => 'FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        'Mitochondrial-Echinoderm' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        'Euplotid' => "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Candida' => "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Ascidian' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        'Mitochondrial-Chlorophycean' => "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Trematode' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        'Mitochondrial-Scenedesmus_obliquus' => "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Thraustochytrium' => "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mitochondrial-Pterobranchia' => "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        'SR1_Gracilibacteria' => "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Pachysolen_tannophilus' => "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Mesodinium' => "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        'Peritrich' => "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        );

    my $exit_code = 0;
    
    foreach my $genetic_code (sort keys %expected_translation_tables) {
        my $expected_translation = $expected_translation_tables{$genetic_code};

        my $reconstructed_translation_table = &show_translation_table($genetic_code);
        my @pts = split(/\n/, $reconstructed_translation_table);
        my $trans_table = shift @pts;

        if ($trans_table eq $expected_translation) {
            print "$genetic_code\tOK\n";
        }
        else {
            print "\n$genetic_code\tERROR:\n"
                . "$expected_translation (EXPECTED)\n"
                . "$trans_table (Computed)\n\n";
            
            $exit_code = 1;
            
        }
        
    }

    exit($exit_code);
    
}

unless (caller) {

    my $genetic_code = $ARGV[0];

    if ($genetic_code eq "TEST") {
        &run_test();
    }

    if ($genetic_code) {
        print "\nTranslation table for $genetic_code:\n\n";
        print &show_translation_table($genetic_code) . "\n\n";
    }
    else {
        print "Supported genetic codes:\n\n" . join("\n", &get_genetic_codes()) . "\n\n"
            . "Choose one to show translation table.\n\n";
    }

    exit(0);
    
}



1; #end of module




