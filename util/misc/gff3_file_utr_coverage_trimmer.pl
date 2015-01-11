#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;
use Fasta_reader;
use GFF3_utils;
use Carp;
use Nuc_translator;
use Data::Dumper;

my $usage = "\n\nusage: $0 gff3_file utr_trim_dat [DEBUG]\n\n";

my $gff3_file = $ARGV[0] or die $usage;
my $utr_trim_dat = $ARGV[1] or die $usage;
my $DEBUG = $ARGV[2] || 0;

my %utr_trim_info;
{
    open (my $fh, $utr_trim_dat) or die $!;
    while (<$fh>) {
        chomp;
        my ($acc, $coords, @rest) = split(/\t/);
        my ($lend, $rend) = split(/-/, $coords);
        $utr_trim_info{$acc} = [$lend, $rend];
    }
    close $fh;
}


my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

        my $trimmed_flag = 0;
		
        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
 
            if ($isoform->has_UTRs()) {
                
                my $isoform_acc = $isoform->{Model_feat_name};
                
                if (my $trim_coords_aref = $utr_trim_info{$isoform_acc}) {
                    
                    my $cdna_len = $isoform->get_cDNA_length();
                    
                    my $targeted_cdna_len = $trim_coords_aref->[1] - $trim_coords_aref->[0] + 1;
                    

                    print "\n====\nBefore:\n" . $isoform->to_GFF3_format() . "\n" if $DEBUG;
                    $isoform = &trim_isoform($trim_coords_aref, $isoform);
                    
                    my $new_cdna_len = $isoform->get_cDNA_length();
                    print "\nLengths, before: $cdna_len, targeted: $targeted_cdna_len, now: $new_cdna_len\n" if $DEBUG;
                    
                    print "\nAfter: (trimmed cdna length: " . join("-", @$trim_coords_aref) . " of len $cdna_len\n" if $DEBUG;
                    print "# $isoform_acc trimmed from $cdna_len to $new_cdna_len cdna length\n";
                    
                    $trimmed_flag = 1;

                }
            }
                        
        }

        
        if ($trimmed_flag) {
            $gene_obj_ref->refine_gene_object();
        }
        
        print $gene_obj_ref->to_GFF3_format() . "\n";
        

    }
}


exit(0);


####
sub trim_isoform {
    my ($trim_coords_aref, $isoform) = @_;

    my ($trim_end5, $trim_end3) = @$trim_coords_aref;
    
    my $cdna_length = 0;

    my $orient = $isoform->get_orientation();
    my @exons = $isoform->get_exons();

    @exons = sort {$a->{end5}<=>$b->{end3}} @exons;

    my @coordsets;

    for (my $i = 0; $i <= $#exons; $i++) {
        
        my $exon = $exons[$i];
        
        my ($end5, $end3) = $exon->get_coords();
        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
        
        $cdna_length += $rend - $lend + 1;
        
        if (my $cds = $exon->get_CDS_obj()) {
            
            my ($cds_end5, $cds_end3) = $cds->get_coords();
            my ($cds_lend, $cds_rend) = sort {$a<=>$b} ($cds_end5, $cds_end3);
            
            push (@coordsets, [ [$lend,$rend], [$cds_lend,$cds_rend] ]);
        }
        else {
            
            push (@coordsets, [ [$lend,$rend], [] ]);
        }
    }
    

    my ($trim_lend, $trim_rend) = ($trim_end5, $trim_end3);

    if ($orient eq '-') {
        $trim_lend = $cdna_length - $trim_lend + 1;
        $trim_rend = $cdna_length - $trim_rend + 1;
        ($trim_lend, $trim_rend) = ($trim_rend, $trim_lend);
    }
    

    if ($DEBUG) {
        print Dumper(\@coordsets);
        print "Trimming coordinates: $trim_lend - $trim_rend\n";
    }
    


    ## do trimming.
    
    my @new_coordsets;
    my $cdna_lend_length = 0;
    foreach my $coordset (@coordsets) {
        my ($exon_coordset, $cds_coordset) = @$coordset;

        my ($exon_lend, $exon_rend) = @$exon_coordset;
        my ($cds_lend, $cds_rend) = @$cds_coordset;
        
        my ($original_exon_lend, $original_exon_rend) = ($exon_lend, $exon_rend);

        my $exon_len = $exon_rend - $exon_lend + 1;
        
        ## checking Left end:
        {
            if ((! $cds_lend) && $exon_len + $cdna_lend_length < $trim_lend) {
                # utr exon trimmed off
                # erase
                ($exon_lend, $exon_rend) = (undef, undef);
            }
            elsif ($trim_lend > $cdna_lend_length && $trim_lend <= $cdna_lend_length + $exon_len) {
                ## trim point within exon
                my $delta = $trim_lend - $cdna_lend_length;
                my $new_exon_lend = $exon_lend + $delta - 1;
                $exon_lend = $new_exon_lend;
                if ($cds_lend && $exon_lend > $cds_lend) {
                    $exon_lend = $cds_lend; # just erasing the left utr
                }
            }
        }
        
        ## checking right end
        if ($exon_lend) { 
            
            if ( (! $cds_lend) && $trim_rend < $cdna_lend_length) {

                # trim point is before utr exon
                # erase it
                ($exon_lend, $exon_rend) = (undef, undef);
            }
            elsif ($trim_rend > $cdna_lend_length && $trim_rend <= $cdna_lend_length + $exon_len) {

                my $delta = $trim_rend - $cdna_lend_length;
                my $new_exon_rend = $original_exon_lend + $delta - 1;
                $exon_rend = $new_exon_rend;
                
                if ($cds_lend && $exon_rend < $cds_rend) {
                    # just removing right utr
                    $exon_rend = $cds_rend;
                }
            }
            
        }
        
        push (@new_coordsets, [ [$exon_lend, $exon_rend], [$cds_lend, $cds_rend] ]);
        
        $cdna_lend_length += $exon_len;
        
    }
    
    ## update gene object coordinates.
    my %exon_coords;
    my %cds_coords;
    foreach my $new_coordset (@new_coordsets) {
        my ($exon_coords_aref, $cds_coords_aref) = @$new_coordset;
        
        my ($exon_lend, $exon_rend) = @$exon_coords_aref;
        my ($cds_lend, $cds_rend) = @$cds_coords_aref;

        if ($exon_lend) {
            my ($exon_end5, $exon_end3) = ($orient eq '+') ? ($exon_lend, $exon_rend) : ($exon_rend, $exon_lend);
            $exon_coords{$exon_end5} = $exon_end3;
        }
        if ($cds_lend) {
            my ($cds_end5, $cds_end3) = ($orient eq '+') ? ($cds_lend, $cds_rend) : ($cds_rend, $cds_lend);
            $cds_coords{$cds_end5} = $cds_end3;
        }

    }

    
    $isoform->populate_gene_obj(\%cds_coords, \%exon_coords);

    return($isoform);
}


                

