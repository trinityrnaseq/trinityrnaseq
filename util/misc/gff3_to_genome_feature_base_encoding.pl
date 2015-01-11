#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;
use GFF3_utils;
use Data::Dumper;

my $usage = "usage: $0 file.gff3\n\n";
my $gff3_file = $ARGV[0] or die $usage;

my %feat_types = ( intergenic => 0,
                   intron => 1,
                   'exon+' => 2,
                   'exon-' => 3,
                   'rRNA+' => 4,
                   'rRNA-' => 5,
    );

my $gene_obj_indexer_href = {};

## associate gene identifiers with contig id's.
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @pos_vec;

    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

        foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
            my @introns = $isoform->get_intron_coordinates();
            
            foreach my $intron (@introns) {
                my ($lend, $rend) = sort {$a<=>$b} @$intron;
                
                for (my $i = $lend; $i <= $rend; $i++) {

                    if ( (! defined $pos_vec[$i]) || $pos_vec[$i] < $feat_types{intron}) {
                        
                        $pos_vec[$i] = $feat_types{intron};
                    }
                }
            }
            
            my $orient = $isoform->get_orientation();
            my $feat_type = "exon$orient";
             if ($isoform->{com_name} =~ /rrna/i && ! $isoform->has_CDS()) {
                $feat_type = "rRNA$orient";
            }
            
            my @exons = $isoform->get_exons();
            foreach my $exon (@exons) {
                
                my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
                for (my $i = $lend; $i <= $rend; $i++) {
                    
                    if ( (! defined $pos_vec[$i]) || $pos_vec[$i] < $feat_types{$feat_type}) {
                        
                        $pos_vec[$i] = $feat_types{$feat_type};
                    }
                }
            }

        }
    }

    shift @pos_vec; #rid first position, since we were using 1-based coordinates above.
    
    foreach my $pos (@pos_vec) {
        unless (defined $pos) {
            $pos = 0;
        }
    }

    my $pos_string = join("", @pos_vec);
    $pos_string =~ s/(\S{60})/$1\n/g;
    
    print ">$asmbl_id\n$pos_string\n";
}


exit(0);

