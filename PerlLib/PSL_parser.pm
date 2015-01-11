package PSL_parser;

use strict;
use warnings;

use Carp;

sub new {
	my $packagename = shift;

	my ($psl_file) = @_;
	
	my $fh;
	
	if (ref $psl_file eq "IO::Handle") {
		$fh = $psl_file; # a filehandle not a file
	}
	else {
	

		unless (-e $psl_file) { 
			confess "error, cannot find $psl_file";
		}
		open ($fh, $psl_file) or confess "Error, cannot open file $psl_file";
	}


	my $self = { file => $psl_file,
				 fh => $fh,
	};

	bless ($self, $packagename);

	
	
	return($self);
}

sub get_next {
	my $self = shift;
	
	my $fh = $self->{fh};
	
	while (my $line = <$fh>) {
		if ($line =~ /^\d+\s/) {
			return(PSL_entry->new($line));
		}
	}

	return(undef); # no more lines
}



##################################################################################
package PSL_entry;

use strict;
use warnings;
use Carp;


sub new {
	my $packagename = shift;
	my ($psl_line) = @_;

	my @fields = split(/\t/, $psl_line);

	my $self = { fields => [@fields] };

	bless ($self, $packagename);

	return($self);
}



################
# blat format: # Q=cDNA T=genomic
################

#  0: match
#  1: mis-match
#  2: rep. match
#  3: N's
#  4: Q gap count
#  5: Q gap bases
#  6: T gap count
#  7: T gap bases
#  8: strand
#  9: Q name
# 10: Q size
# 11: Q start
# 12: Q end
# 13: T name
# 14: T size
# 15: T start
# 16: T end
# 17: block count
# 18: block Sizes
# 19: Q starts
# 20: T starts
# 21: Q seqs (pslx format)
# 22: T seqs (pslx format)

## All sequences start at 0 here; array-based.


sub get_line {
	my $self = shift;
	return(join("\t", @{$self->{fields}}));
}


sub get_match_count {
	my $self = shift;
	return($self->{fields}->[0]);
}

sub get_mismatch_count {
	my $self = shift;
	return($self->{fields}->[1]);
}

sub get_N_count {
	my $self = shift;
	return($self->{fields}->[3]);
}

sub get_Q_gap_count {
	my $self = shift;
	return($self->{fields}->[4]);
}

sub get_Q_gap_bases {
	my $self = shift;
	return($self->{fields}->[5]);
}

sub get_T_gap_count {
	my $self = shift;
	return($self->{fields}->[6]);
}

sub get_T_gap_bases {
	my $self = shift;
	return($self->{fields}->[7]);
}

sub get_strand {
	my $self = shift;
	return($self->{fields}->[8]);
}

sub get_Q_name {
	my $self = shift;
	return($self->{fields}->[9]);
}

sub get_Q_size {
	my $self = shift;
	return($self->{fields}->[10]);
}

sub get_Q_span {
	my $self = shift;
	return($self->{fields}->[11] + 1, $self->{fields}->[12]);
}

sub get_T_name {
	my $self = shift;
	return($self->{fields}->[13]);
}

sub get_T_size {
	my $self = shift;
	return($self->{fields}->[14]);
}

sub get_T_span {
	my $self = shift;
	return($self->{fields}->[15] + 1, $self->{fields}->[16]);
}


####
sub get_per_id {
	my $self = shift;

	my $matches = $self->get_match_count();
	my $mismatches = $self->get_mismatch_count();

	my $per_id = $matches / ($matches + $mismatches) * 100;

	$per_id = sprintf("%.2f", $per_id);
	
	return($per_id);
}




####
sub get_alignment_coords {
	my $self = shift;

	my @x = @{$self->{fields}};
	

    my @alignment_segments;
 
 
    my @cdna_coords = split (/,/, $x[19]);
    my @genomic_coords = split (/,/, $x[20]);
    my @lengths = split (/,/, $x[18]);
 
	my $strand = $self->get_strand();
	my $cdna_length = $self->get_Q_size();
	
	my @ret_genome_coords;
	my @ret_cdna_coords;

    ## report each segment match as a separate btab entry:
    my $segment_number = 0;
    my $num_segs = scalar(@genomic_coords);
	for (my $i = 0; $i < $num_segs; $i++) {
		$segment_number++;
		my $length = $lengths[$i];
		unless (defined $length) { next; }
		my $cdna_coord = $cdna_coords[$i];
		my $genomic_coord = $genomic_coords[$i];
		my ($cdna_end5, $cdna_end3) = (++$cdna_coord, $cdna_coord + $length - 1);
		my ($genomic_end5, $genomic_end3) = (++$genomic_coord, $genomic_coord + $length -1);
		if ($strand eq "-") {
			($cdna_end5, $cdna_end3) = ($cdna_length - $cdna_end5 + 1, $cdna_length - $cdna_end3 + 1);
		}
		
		push (@ret_genome_coords, [$genomic_end5, $genomic_end3]);
		push (@ret_cdna_coords, [$cdna_end5, $cdna_end3]);
	}
	

	return(\@ret_genome_coords, \@ret_cdna_coords);
}


sub toString {
	my $self = shift;

	my $genome_acc = $self->get_T_name();
	my $cdna_acc = $self->get_Q_name();

	my $genome_length = $self->get_T_size();
	my $cdna_length = $self->get_Q_size();
	

	my $strand = $self->get_strand();

	my ($genome_lend, $genome_rend) = $self->get_T_span();
	my ($cdna_lend, $cdna_rend) = $self->get_Q_span();
	
	my $per_id = sprintf("%.2f", $self->get_per_id());

	my $ret_text = join("\t", $genome_acc, "($genome_length)", "$genome_lend-$genome_rend", $strand,
						$cdna_acc, "($cdna_length)", "$cdna_lend-$cdna_rend", $per_id);
	

	my ($genome_coords_aref, $cdna_coords_aref) = $self->get_alignment_coords();
	
	my @genome_coords = @$genome_coords_aref;
	my @cdna_coords = @$cdna_coords_aref;

	my $align_text = "";
	while (@genome_coords) {
		my $genome_coordset = shift @genome_coords;
		my $cdna_coordset = shift @cdna_coords;
		
		if ($align_text) {
			$align_text .= "....";
		}
		$align_text .= $genome_coordset->[0] . "(" . $cdna_coordset->[0] . ")-"
			. $genome_coordset->[1] . "(" . $cdna_coordset->[1] . ")";
	}
	
	$ret_text .= "\n$align_text\n";
	
	return ($ret_text);
}


1; #EOM






