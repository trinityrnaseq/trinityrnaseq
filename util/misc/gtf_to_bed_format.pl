#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $usage = "usage: $0 transcripts.gtf\n\n";

my $gtf_file = $ARGV[0] or die $usage;


main: {
	
	my %genome_trans_to_coords;
	
	open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
	while (<$fh>) {
		chomp;
		
		unless (/\w/) { next; }
		
		my @x = split(/\t/);
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];

		my $orient = $x[6];
		
		my $info = $x[8];
		
		unless ($type eq 'exon') { next; }

		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/\s+/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{gene_id} or die "Error, no gene_id at $_";
		my $trans_id = $atts{transcript_id} or die "Error, no trans_id at $_";
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$genome_trans_to_coords{$scaff}->{$gene_id}->{$trans_id}->{$end5} = $end3;

	}


	## Output genes in gff3 format:

	foreach my $scaff (sort keys %genome_trans_to_coords) {

		my $genes_href = $genome_trans_to_coords{$scaff};

		foreach my $gene_id (keys %$genes_href) {

			my $trans_href = $genes_href->{$gene_id};

			foreach my $trans_id (keys %$trans_href) {

				my $coords_href = $trans_href->{$trans_id};

				my $gene_obj = new Gene_obj();

				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = $trans_id;
				$gene_obj->{com_name} = "$gene_id $trans_id";
				
				$gene_obj->{asmbl_id} = $scaff;
				
				$gene_obj->populate_gene_object($coords_href, $coords_href);
			
				print $gene_obj->to_BED_format();
								
			}
		}
	}


	exit(0);
}

