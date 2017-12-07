#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;
use Nuc_translator;
use Data::Dumper;

use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

###################################################
#
# Required:
#
#  --genome_fa <string>          genome fasta file
# 
#  --flattened_gff <string>      flattened gff file
#
# Optional:
#
#  --no_revcomp                  do not revcomp reverse strand transcripts
# 
###################################################


__EOUSAGE__

    ;


my $genome_fa;
my $flattened_gff;
my $NO_REVCOMP_FLAG = 0;
my $help_flag;

&GetOptions ( 'h' => \$help_flag,

              'genome_fa=s' => \$genome_fa,
              'flattened_gff=s' => \$flattened_gff,

              'no_revcomp' => \$NO_REVCOMP_FLAG,
    );


if ($help_flag) {
    die $usage;
}

unless ($genome_fa && $flattened_gff) {
    die $usage;
}

main: {

    
    my $fasta_reader = new Fasta_reader($genome_fa);
    my %seqs = $fasta_reader->retrieve_all_seqs_hash();

    my %gene_to_transcripts = &parse_gff($flattened_gff);

    unless (%gene_to_transcripts) {
        die "Error, no gff records parsed.  Be sure to use the flattened gff file";
    }
    
    #die Dumper(\%gene_to_transcripts);
    
    foreach my $gene (keys %gene_to_transcripts) {

        my @transcripts = keys %{$gene_to_transcripts{$gene}};
        my $transcript_counter = 0;
        foreach my $transcript (@transcripts) {
            $transcript_counter += 1;
            my @regions = @{$gene_to_transcripts{$gene}->{$transcript}};
            @regions = sort {$a->{lend}<=>$b->{lend}} @regions;
            my $trans_seq = "";
            my @path_coords;
            foreach my $region (@regions) {
                my $lend = $region->{lend};
                my $rend = $region->{rend};
                my $chr = $region->{chr};
                my $exon_number = $region->{exon_number};
                my $seg_len = $rend - $lend + 1;
                my $seq_seg = substr($seqs{$chr}, $lend-1, $seg_len);
                
                push (@path_coords, [$exon_number, length($trans_seq)+1, length($trans_seq) + length($seq_seg)]);
                $trans_seq .= $seq_seg;
            }
            if ($regions[0]->{orient} eq '-' && ! $NO_REVCOMP_FLAG) {
                # revcomp everthing
                $trans_seq = &reverse_complement($trans_seq);
                my $trans_len = length($trans_seq);
                foreach my $path_coord (@path_coords) {
                    my $new_rend = $trans_len - $path_coord->[1] + 1;
                    my $new_lend = $trans_len - $path_coord->[2] + 1;

                    ($path_coord->[1], $path_coord->[2]) = ($new_lend, $new_rend);
                }
                @path_coords = sort {$a->[1]<=>$b->[1]} @path_coords;
            }

            ## construct header
            my @path_text;
            foreach my $path_coord (@path_coords) {
                my ($exon_number, $lend, $rend) = @$path_coord;
                $lend--;
                $rend--;
                push (@path_text, "$exon_number:$lend-$rend");
            }
            
            
            my $header = ">${gene}_i${transcript_counter} path=[" . join(" ", @path_text) . "]";

            print "$header\n$trans_seq\n";
        }
    }

    exit(0);
    
              
}


####
sub parse_gff {
    my ($gff_file) = @_;

    my %gene_to_transcripts;
    
    open(my $fh, $gff_file) or die $!;
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my $chr = $x[0];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $orient = $x[6];
        my $info = $x[8];

        unless ($feat_type eq "exonic_part") { next; }
        
        $info =~ /transcripts \"([^\"]+)\";.*\s+exonic_part_number \"(\d+)\"; gene_id \"([^\"]+)\"/ or die "Error, cannot parse $info";

        my $transcripts_list = $1;
        my $exonic_part_number = $2;
        my $gene_id = $3;

        my $region_struct = { chr => $chr,
                              lend  => $lend,
                              rend => $rend,
                              orient => $orient,
                              exon_number => $exonic_part_number,
        };

        my @transcripts = split(/\+/, $transcripts_list);
        foreach my $transcript (@transcripts) {
            push (@{$gene_to_transcripts{$gene_id}->{$transcript}}, $region_struct);
        }
    }

    close $fh;

    return(%gene_to_transcripts);
}
        
        
