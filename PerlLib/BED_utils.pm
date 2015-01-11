package BED_utils;

use strict;
use warnings;
use Carp;
use Gene_obj;

sub index_BED_as_gene_objs {
	my ($gff_filename, $gene_id_to_gene_obj_href) = @_;

	my %contig_to_gene_list;

	open (my $fh, $gff_filename) or die "Error, cannot open file $gff_filename";
	while (<$fh>) {
		if (/^\#/) { next; }
        chomp;
		unless (/\w/) { next; }
		
		my $bed_line = $_;
		
		my $gene_obj;

		eval {
			$gene_obj = &Gene_obj::BED_line_to_gene_obj($bed_line);
			my @introns = $gene_obj->get_intron_coordinates(); # this method breaks if all exons are single bases.  Ignore these weird things.
		};

		if ($@) {
			print STDERR "ERROR, cannot create gene for bed line:\n$bed_line\n$@\n";
			next;
		}
		

		my $gene_id = $gene_obj->{TU_feat_name};
		
		my $indexed_gene_obj = $gene_id_to_gene_obj_href->{$gene_id};
	    if ($indexed_gene_obj) {
			$indexed_gene_obj->add_isoform($gene_obj);
		}
		else {
			$gene_id_to_gene_obj_href->{$gene_id} = $gene_obj;
			my $contig = $gene_obj->{asmbl_id};
			push (@{$contig_to_gene_list{$contig}}, $gene_id);
		}
	}
	close $fh;

	return(\%contig_to_gene_list);
}



1; #EOM
	
