#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 fileA.maps fileB.maps gene_annots.gff3\n\n";

my $file_A_maps = $ARGV[0] or die $usage;
my $file_B_maps = $ARGV[1] or die $usage;
my $gene_annots_gff3 = $ARGV[2] or die $usage;

main: {

    my ($A_fusion_to_genes_href, $A_FL_genes) = &parse_map_file($file_A_maps);
    
    my ($B_fusion_to_genes_href, $B_FL_genes) = &parse_map_file($file_B_maps);

    my %gene_structs = &get_gene_positions($gene_annots_gff3);
        
    my %all_fusions = map { $_ => 1 } (keys %$A_fusion_to_genes_href, keys %$B_fusion_to_genes_href);

    foreach my $fusion (keys %all_fusions) {

        my $neighboring = &are_neighboring($fusion, \%gene_structs);

        print "$fusion\t$neighboring\t";
        
        if ($A_fusion_to_genes_href->{$fusion} && $B_fusion_to_genes_href->{$fusion}) {
            print "BOTH\n";
        }
        elsif ($A_fusion_to_genes_href->{$fusion}) {
            print "A_only\n";
        }
        elsif ($B_fusion_to_genes_href->{$fusion}) {
            print "B_only\n";
        }
        
    }
    
    exit(0);
    
}

####
sub parse_map_file {
    my ($map_file) = @_;

    my %fusion_to_genes;
    my %FL_genes;

    open (my $fh, $map_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($fl_entries, $rest) = split(/\t/);
        
        my @genes = split(/,/, $fl_entries);
        
        foreach my $gene (@genes) {

            $FL_genes{$gene} = $fl_entries;
        }

        if (scalar @genes > 1) {
            $fusion_to_genes{$fl_entries} = \@genes;
        }
    }
    close $fh;

    return(\%fusion_to_genes, \%FL_genes);
}


####
sub get_gene_positions {
    my ($annot_file) = @_;

    my %gene_structs;
    my %scaff_to_structs;

    open (my $fh, $annot_file) or die $!;
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        my @x = split(/\t/);
        if ($x[2] eq 'gene') {
            $x[8] =~ /ID=([^;]+)/ or die "Error, cannot parse ID from $x[8]";
            my $gene_id = $1;
            my $struct = { scaffold => $x[0],
                           midpt =>  ($x[3]+$x[4])/2,
                           order => undef,
                           gene_id => $gene_id,
                       };
        
            $gene_structs{$gene_id} = $struct;
            push (@{$scaff_to_structs{$x[0]}}, $struct);
            
        }
    }
    close $fh;


    ## order the genes.
    
    my $pos = 0;
    foreach my $scaff (keys %scaff_to_structs) {
        my @structs = @{$scaff_to_structs{$scaff}};
        @structs = sort {$a->{midpt}<=>$b->{midpt}} @structs;

        foreach my $struct (@structs) {
            $pos++;
            $struct->{order} = $pos;
        }
        $pos++; # gap between scaffs
    }

    return(%gene_structs);
}

####
sub are_neighboring {
    my ($fusion, $gene_structs_href) = @_;

    my @gene_structs;
    foreach my $gene (split(/,/, $fusion)) {
        my ($trans_id, $gene_id) = split(/;/, $gene);
        
        my $gene_struct = $gene_structs_href->{$gene_id} or die "Error, no gene struct for $gene_id";
        
        push (@gene_structs, $gene_struct);
    }
    
    @gene_structs = sort {$a->{order}<=>$b->{order}} @gene_structs;
    
    if ($gene_structs[1]->{order} - $gene_structs[0]->{order} <= 4 ) {
        return("Neighbors");
    }
    else {
        return("NOT_neighbors ($gene_structs[0]->{order} vs. $gene_structs[1]->{order})");
    }
}
