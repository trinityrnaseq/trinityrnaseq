package SAM_to_AlignGraph;

## Static class

use strict;
use warnings;

use AlignGraph; 
use Carp;

use SAM_reader;
use SAM_entry;

sub construct_AlignGraph {
	my ($sam_file) = @_;

	my $graph = new AlignGraph();

	my $sam_reader = new SAM_reader($sam_file);
	
	my $counter = 0;

	while ($sam_reader->has_next()) {

		$counter++;
		print STDERR "\r[$counter]   " if $counter % 100 == 0;

		my $sam_entry = $sam_reader->get_next();

		my $scaff = $sam_entry->get_scaffold_name();
		my $read_acc = $sam_entry->get_read_name();
		
		if ($sam_entry->is_query_unmapped()) { next; }
		
		
		
		my $query_strand = $sam_entry->get_query_transcribed_strand();

		my ($genome_coords_aref, $query_coords_aref) = $sam_entry->get_alignment_coords();
				
		$graph->add_alignment($read_acc, $scaff, $query_strand, $genome_coords_aref);
		
	}

	return($graph);
}

1; #EOM


