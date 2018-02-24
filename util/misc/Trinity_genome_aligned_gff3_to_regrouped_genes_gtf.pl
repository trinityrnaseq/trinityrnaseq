#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Data::Dumper;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use SingleLinkageClusterer;
use Overlap_piler;

my $usage = "usage: $0 gmap.gff3 Trinotate.annot_mappings\n\n";

my $gmap_gff3 = $ARGV[0] or die $usage;
my $trinotate_mapping_file = $ARGV[1] or die $usage;


my $DEBUG = 0;

my $GENERIC_GENE_COUNTER = 0;

my %GENE_NAMES_USED;

main: {

    my %annotations = &parse_annotation_mappings($trinotate_mapping_file);
    
    my %contig_to_transcripts = &parse_gmap_gff3($gmap_gff3);

    foreach my $contig (keys %contig_to_transcripts) {
        
         my @transcript_structs = values %{$contig_to_transcripts{$contig}};

         # order the exon coords
         foreach my $struct (@transcript_structs) {
             @{$struct->{coords}} = sort {$a->[0]<=>$b->[0]} @{$struct->{coords}};
         }
         
         my ($plus_structs_aref, $minus_structs_aref) = &separate_by_strand(@transcript_structs);

         foreach my $struct_list_aref ($plus_structs_aref, $minus_structs_aref) {

             unless (@$struct_list_aref) { next; } # nothing to do
             
             my @gene_grouped_clusters = &group_genes_by_span_overlaps($struct_list_aref);

             @gene_grouped_clusters = &group_genes_by_exon_overlaps(@gene_grouped_clusters);
                          
             
             foreach my $gene_grouped_cluster (@gene_grouped_clusters) {
                 print STDERR "Contig: $contig, " . Dumper($gene_grouped_cluster) if $DEBUG;

                 my ($gene_id_orig, $gene_id_use)  = &get_gene_id($gene_grouped_cluster);

                 my $gene_name = $annotations{$gene_id_orig} || "$gene_id_orig";
                 
                 foreach my $struct (@$gene_grouped_cluster) {
                     my $strand = $struct->{strand};
                     my $align_name = $struct->{align_name};
                     my $transcript_name = $align_name;
                     $transcript_name =~ s/.mrna\d+//;

                     my $trans_annot = $annotations{$transcript_name} || "$transcript_name";
                     
                     my @coords = @{$struct->{coords}};
                     foreach my $coordset (@coords) {
                         print join("\t", $contig, ".", "exon", $coordset->[0], $coordset->[1], ".", $strand, ".",
                                    "transcript_id \"$align_name\"; gene_id \"$gene_id_use\"; transcript_biotype \"processed_transcript\"; gene_name \"$gene_name\"; transcript_name \"$transcript_name\"") . "\n";
                     }
                     print "\n"; # spacer
                 }
                 
             }

         }
         
         
     }
 

     exit(0);
     
}

####
sub parse_gmap_gff3 {
    my ($gmap_gff3_file) = @_;

    ##  Example record
    ##
    # AMEXG_0030000816        AmexG_v3.0.0.fa.gmap    gene    4572879 4573094 .       -       .       ID=c104_g1_i1.path1;Name=c104_g1_i1
    # AMEXG_0030000816        AmexG_v3.0.0.fa.gmap    mRNA    4572879 4573094 .       -       .       ID=c104_g1_i1.mrna1;Name=c104_g1_i1;Parent=c104_g1_i1.path1;coverage=100.0;identity=100.0;matches=216;mismatches=0;indels=0;unknowns=0
    # AMEXG_0030000816        AmexG_v3.0.0.fa.gmap    exon    4572879 4573094 100     -       .       ID=c104_g1_i1.mrna1.exon1;Name=c104_g1_i1;Parent=c104_g1_i1.mrna1;Target=c104_g1_i1 1 216 +
    # AMEXG_0030000816        AmexG_v3.0.0.fa.gmap    CDS     4572880 4573092 100     -       0       ID=c104_g1_i1.mrna1.cds1;Name=c104_g1_i1;Parent=c104_g1_i1.mrna1;Target=c104_g1_i1 3 215 +
        

    my %contig_to_transcripts;
    
    ## just capture the exon records.
    open(my $fh, $gmap_gff3_file) or die $!;
    while(<$fh>) {
        if (/^\#/) { next; }
        unless (/\w/) { next; }
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my $contig_id = $x[0];
        my $feat_type = $x[2];
        my $lend = $x[3];
        my $rend = $x[4];
        my $strand = $x[6];
        my $info = $x[8];

        
        unless ($feat_type eq "exon") { next; }

        my $align_name;
        if ($info =~ /Parent=([^;]+)/) {
            $align_name = $1;
        }
        else{
            die "Error, cannot extract alignment name (Parent) from info: $info of line: $line";
        }

        push (@{$contig_to_transcripts{$contig_id}->{$align_name}->{coords}}, [$lend, $rend]);
        
        $contig_to_transcripts{$contig_id}->{$align_name}->{strand} = $strand;
        $contig_to_transcripts{$contig_id}->{$align_name}->{align_name} = $align_name;
        
    }
    close $fh;
 

    return(%contig_to_transcripts);
    
}

####
sub separate_by_strand {
    my (@transcript_structs) = @_;

    my @plus_structs;
    my @minus_structs;
    
    foreach my $struct (@transcript_structs) {
        if ($struct->{strand} eq '+') {
            push (@plus_structs, $struct);
        }
        elsif ($struct->{strand} eq '-') {
            push (@minus_structs, $struct);
        }
        else {
            die "Error, not sure what strand this corresponds to: " . Dumper($struct);
        }
    }

    return(\@plus_structs, \@minus_structs);
}

####
sub group_genes_by_span_overlaps {
    my ($struct_list_aref) = @_;

    print STDERR "\n\n// Grouping by overlaps: " . Dumper($struct_list_aref) . "\n" if $DEBUG;
    
    my %id_to_struct;

    my $piler = new Overlap_piler();
    
    foreach my $struct (@$struct_list_aref) {
        my $align_name = $struct->{align_name};

        $id_to_struct{$align_name} = $struct;

        my @coordsets = @{$struct->{coords}};

        my @all_coords;
        foreach my $coordset (@coordsets) {
            my ($lend, $rend) = @$coordset;
            push (@all_coords, $lend, $rend);
        }
        @all_coords = sort {$a<=>$b} @all_coords;

        my $span_lend = shift @all_coords;
        my $span_rend = pop @all_coords;

        
        $piler->add_coordSet($align_name, $span_lend, $span_rend);
        print STDERR "-adding $align_name, $span_lend-$span_rend\n" if $DEBUG;
    }

    my @clusters = $piler->build_clusters();

    print STDERR "Clusters: " . Dumper(\@clusters) . "\n" if $DEBUG;

    my @struct_clusters;
    foreach my $cluster (@clusters) {
        my @ids = @$cluster;

        my @struct_cluster;
        foreach my $id (@ids) {
            my $struct = $id_to_struct{$id};
            push (@struct_cluster, $struct);
        }
        push (@struct_clusters, [@struct_cluster]);
    }

    return(@struct_clusters);

}

####
sub get_gene_id {
    my ($struct_list_aref) = @_;

    my %genes;

    foreach my $struct (@$struct_list_aref) {

        my $align_name = $struct->{align_name};
        $align_name =~ s/_i\d+.mrna\d+//;

        $genes{$align_name}++;
    }

    my @gene_ids = keys %genes;
    if (scalar(@gene_ids) == 1 && ! exists $GENE_NAMES_USED{$gene_ids[0]}) {
        return($gene_ids[0], $gene_ids[0]);
    }
    else {
        $GENERIC_GENE_COUNTER++;
        return($gene_ids[0], "G_$GENERIC_GENE_COUNTER");
    }
}


####
sub parse_annotation_mappings {
    my ($trinotate_mapping_file) = @_;

    my %annots;
    
    open(my $fh, $trinotate_mapping_file) or die $!;
    while(<$fh>) {
        chomp;
        my ($feature_id, $annot) = split(/\t/);
        $annots{$feature_id} = $annot;
    }
    close $fh;

    return(%annots);
}

####
sub group_genes_by_exon_overlaps {
    my (@groups) = @_;

    my @groups_ret = ();

    foreach my $group (@groups) {
        if (scalar(@$group) == 1) {
            push (@groups_ret, $group);
        }
        else {
            my @refined_groups = &partition_by_exon_overlap(@$group);
            push (@groups_ret, @refined_groups);
        }
    }

    return(@groups_ret);
}

####
sub partition_by_exon_overlap {
    my @structs = @_;

    my %id_to_struct;
    foreach my $struct (@structs) {
        $id_to_struct{$struct->{align_name}} = $struct;
    }
    
    my @pairs;

    for(my $i = 0; $i < $#structs; $i++) {
        my $struct_i = $structs[$i];
        
        for (my $j = $i + 1; $j <= $#structs; $j++) {
            my $struct_j = $structs[$j];
            
            if (&exon_overlap($struct_i, $struct_j)) {
                
                push (@pairs, [$struct_i->{align_name}, $struct_j->{align_name}]);
            }
        }
    }

    
    my @clusters;
    if (scalar(@pairs) > 1) {
        @clusters =  &SingleLinkageClusterer::build_clusters(@pairs);
    }
    elsif (scalar @clusters == 1) {
        @clusters = @pairs;
    }
    
    my %seen;

    my @clusters_ret;
    foreach my $cluster (@clusters) {
        my @cluster_ret;
        foreach my $ele (@$cluster) {
            my $struct = $id_to_struct{$ele} or die "Error, no struct for $ele";
            push (@cluster_ret, $struct);
            $seen{$ele} = 1;
        }
        push (@clusters_ret, \@cluster_ret);
    }

    foreach my $struct (@structs) {
        if (! exists $seen{$struct->{align_name}}) {
            push (@clusters_ret, [$struct]);
        }
    }

    return(@clusters_ret);
}

####
sub exon_overlap {
    my ($struct_A, $struct_B) = @_;

    print "A: " . Dumper($struct_A) . "\nB: " . Dumper($struct_B) . "\n" if $DEBUG;
    
    my @exon_coords_A = @{$struct_A->{coords}};
    my @exon_coords_B = @{$struct_B->{coords}};

    foreach my $exon_A (@exon_coords_A) {
        foreach my $exon_B (@exon_coords_B) {

            if ($exon_A->[0] < $exon_B->[1] && $exon_A->[1] > $exon_B->[0]) {
                return(1); # overlap detected
            }
        }
    }

    return(0); # no overlap
}

