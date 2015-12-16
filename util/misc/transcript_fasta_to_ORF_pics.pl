#!/usr/bin/env perl

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");

use strict;
use warnings;
use Getopt::Std;
use Fasta_reader;
use Longest_orf;
use List::Util qw (min max);

use Bio::Graphics;
use Bio::SeqFeature::Generic;

my $usage = "\n\nusage: $0 transcripts.fasta\n\n";

my $transcripts_fasta_file = $ARGV[0] or die $usage;
my $min_prot_len = 50;


main: {

    my $fasta_reader = new Fasta_reader($transcripts_fasta_file);
    

    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        
        my $seq_length = length($sequence);
        

        my $panel = Bio::Graphics::Panel->new(
                                              -length    => $seq_length,
                                              -width     => 800,
                                              -pad_left  => 10,
                                              -pad_right => 10,
                                              );
        
        my $full_length = Bio::SeqFeature::Generic->new(
                                                        -start => 1,
                                                        -end   => $seq_length,
                                                        -strand => 1,
                                                        -display_name => $acc,
                                                        );
        
        
        
        ## add ticker
        $panel->add_track($full_length,
                          -glyph   => 'arrow',
                          -tick    => 2,
                          -fgcolor => 'black',
                          -double  => 1,
                          
                          );
        
        ## add feature
        $panel->add_track($full_length,
                          -glyph   => 'transcript2',
                          -bgcolor => 'black',
                          -label => 1,
                          );


        my $longest_orf_finder = new Longest_orf();
        $longest_orf_finder->allow_5prime_partials();
        $longest_orf_finder->allow_3prime_partials();

        my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);
        
        

        
        @orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;

        
        

        my $top_track = $panel->add_track(-glyph => 'transcript2',
                                          -label => 1,
                                          -strand_arrow => 1,
                                          -bgcolor => 'blue',
                                          
                                          );
        
        my $bottom_track = $panel->add_track(-glyph => 'transcript2',
                                             -label => 1,
                                             -strand_arrow => 1,
                                             -bgcolor => 'red',
                                             );
        

        my $pep_file = "$acc.orfs.pep";
        my $cds_file = "$acc.orfs.cds";
        
        open (my $pep_ofh, ">$pep_file") or die "Error, cannot write to $pep_file";
        open (my $cds_ofh, ">$cds_file") or die "Error, cannot write to $cds_file";
        
        foreach my $orf (@orf_structs) {
            
            my $start = $orf->{start};
            my $stop = $orf->{stop};
            
            if ($stop <= 0) { $stop += 3; } # edge issue
            
            my $length = int($orf->{length}/3);
            if ($length < $min_prot_len) { next; }
            
            my $orient = $orf->{orient};
            my $protein = $orf->{protein};
            my $cds = $orf->{sequence};
            

            my $frame = $start % 3;
            if ($frame == 0) {
                $frame = 3;
            }
            
            my $orf_name = "$acc.$start-$stop.F$frame";
            print $pep_ofh ">$orf_name\n$protein\n";
            print $cds_ofh ">$orf_name\n$cds\n";
            

            if ($orient eq '+') {
                
                                
                my $feature = Bio::SeqFeature::Generic->new(
                                                            -start        => $start,
                                                            -end          => $stop,
                                                            -display_name => "orf: $start-$stop (F:$frame)",
                                                            -strand => 1,
                                                            
                                                            );
                $top_track->add_feature($feature);
                
            }
            else {
                # minus strand
            
                $frame = -1 * $frame;
    
                my $feature = Bio::SeqFeature::Generic->new(
                                                            -start        => $start,
                                                            -end          => $stop,
                                                            -display_name => "orf: $start-$stop (F:$frame)",
                                                            -strand => -1,
                                                            
                                                            );
                $bottom_track->add_feature($feature);
            }
        }
        close $pep_ofh;
        close $cds_ofh;
        

        my $img_filename = "$acc.orf.png";
        print STDERR "-writing: $img_filename\n";
        open (my $ofh, ">$img_filename") or die "Error, cannot write to file $img_filename";
        binmode($ofh);
        print $ofh $panel->png();
        close $ofh;
    }
    

    exit(0);
}


