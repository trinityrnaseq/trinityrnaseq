package SAM_entry;

use strict;
use warnings;
use Carp;



sub new {
	my $packagename = shift;
	my ($line) = @_;

	unless (defined $line) {
		confess "Error, need sam text line as parameter";
	}

	chomp $line;

	my @fields = split(/\t/, $line);
	
	my $self = {
		_fields => [@fields],
	};
	
	bless ($self, $packagename);
	
	return($self);
}


####
sub get_read_name {
	my $self = shift;
	return ($self->{_fields}->[0]);
}


####
sub get_scaffold_name {
	my $self = shift;
	return($self->{_fields}->[2]);
}

####
sub get_aligned_position {
	my $self = shift;
	return($self->{_fields}->[3]);
}

sub get_scaffold_position { # preferred
	my $self = shift;
	return($self->get_aligned_position());
}


####
sub get_cigar_alignment {
	my $self = shift;
	return($self->{_fields}->[5]);
}

####
sub get_alignment_coords {
	my $self = shift;

	my $genome_lend = $self->get_aligned_position();

	my $alignment = $self->get_cigar_alignment();

	my $query_lend = 0;

	my @genome_coords;
	my @query_coords;


	$genome_lend--; # move pointer just before first position.
	
	while ($alignment =~ /(\d+)([A-Z])/g) {
		my $len = $1;
		my $code = $2;
		
		unless ($code =~ /^[MSDNIH]$/) {
			confess "Error, cannot parse cigar code [$code]";
		}
		
		# print "parsed $len,$code\n";
		
		if ($code eq 'M' || $code eq 'S' || $code eq 'H') { # aligned bases match or mismatch
			
			my $genome_rend = $genome_lend + $len;
			my $query_rend = $query_lend + $len;
			
			push (@genome_coords, [$genome_lend+1, $genome_rend]);
			push (@query_coords, [$query_lend+1, $query_rend]);
			
			# reset coord pointers
			$genome_lend = $genome_rend;
			$query_lend = $query_rend;
			
		}
		elsif ($code eq 'D' || $code eq 'N') { # insertion in the genome, gap in query (intron, perhaps)
			$genome_lend += $len;
			
		}

		elsif ($code eq 'I') { # gap in genome, insertion in query 

			$query_lend += $len;

		}
	}

	return(\@genome_coords, \@query_coords);
}


####
sub get_mate_scaffold_name {
	my $self = shift;
   
	return($self->{_fields}->[6]);
}


####
sub set_mate_scaffold_name {
	my $self = shift;
	my $mate_scaffold_name = shift;
	
	$self->{_fields}->[6] = $mate_scaffold_name;
	
	return;
}


####
sub get_mate_scaffold_position {
	my $self = shift;

	return($self->{_fields}->[7]);
}


####
sub set_mate_scaffold_position {
	my $self = shift;
	my $scaff_pos = shift;
   
	$self->{_fields}->[7] = $scaff_pos;
	
	return;
}


####
sub toString {
	my $self = shift;
	return( join("\t", @{$self->{_fields}}) );
}


####
sub get_mapping_quality {
	my $self = shift;
	return($self->{_fields}->[4]);
}


####
sub get_sequence {
	my $self = shift;
	return($self->{_fields}->[9]);
}

####
sub get_quality_scores {
	my $self = shift;
	return($self->{_fields}->[10]);
}



###################
## Flag Processing
###################

# from sam format spec:

=flag_description

Flag Description
0x0001 the read is paired in sequencing, no matter whether it is mapped in a pair
0x0002 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
0x0004 the query sequence itself is unmapped
0x0008 the mate is unmapped 1
0x0010 strand of the query (0 for forward; 1 for reverse strand)
0x0020 strand of the mate 1
0x0040 the read is the first read in a pair 1,2
0x0080 the read is the second read in a pair 1,2
0x0100 the alignment is not primary (a read having split hits may have multiple primary alignment records)
0x0200 the read fails platform/vendor quality checks
0x0400 the read is either a PCR duplicate or an optical duplicate

1. Flag 0x02, 0x08, 0x20, 0x40 and 0x80 are only meaningful when flag 0x01 is present.
2. If in a read pair the information on which read is the first in the pair is lost in the upstream analysis, flag 0x01 shuld
be present and 0x40 and 0x80 are both zero.

=cut


####
sub get_flag {
	my $self = shift;
	my $flag = $self->{_fields}->[1];
	return($flag);
}

sub set_flag {
	my $self = shift;
	my $flag = shift;

	unless (defined $flag) {
		confess "Error, need flag value";
	}

	$self->{_fields}->[1] = $flag;
	return;
}

####
sub is_paired {
	my $self = shift;
	return($self->_get_bit_val(0x0001));
}

sub set_paired {
	my $self = shift;
	my $bit_val = shift;
	
	$self->_set_bit_val(0x0001, $bit_val);
	
	return;
}

####
sub is_proper_pair {
	my $self = shift;
	return($self->_get_bit_val(0x0002));
}

sub set_proper_pair {
	my $self = shift;
	my $bit_val = shift;
	
	$self->_set_bit_val(0x0002, $bit_val);
	return;
}

####
sub is_query_unmapped {
	my $self = shift;
	return($self->_get_bit_val(0x0004));
}

sub set_query_unmapped {
	my $self = shift;
	my $bit_val = shift;

	$self->_set_bit_val(0x0004, $bit_val);
}


####
sub is_mate_unmapped {
	my $self = shift;
	return($self->_get_bit_val(0x0008));
}

sub set_mate_unmapped {
	my $self = shift;
	my $bit_val = shift;
	
	return($self->_set_bit_val(0x0008, $bit_val));
}

####
sub get_query_strand {
	my $self = shift;
	
	my $strand = ($self->_get_bit_val(0x0010)) ? '-' : '+';
	return($strand);
}

####
sub get_query_transcribed_strand {
	my $self = shift;

	my $strand = $self->get_query_strand();
	
	if ($self->is_paired() && $self->is_first_in_pair()) {
		
		my $transcribed_strand = ($strand eq '+') ? '-' : '+';
		
		return($transcribed_strand);
	}
	else {
		return($strand);
	}
}




sub set_query_strand {
	my $self = shift;
	my $strand = shift;

	unless ($strand eq '+' || $strand eq '-') {
		confess "Error, strand value must be [+-]";
	}

	my $bit_val = ($strand eq '+') ? 0 : 1;
	$self->_set_bit_val(0x0010, $bit_val);
}

####
sub get_mate_strand {
	my $self = shift;
	
	my $strand = ($self->_get_bit_val(0x0020)) ? '-' : '+';
	return($strand);
}

sub set_mate_strand {
	my $self = shift;
	my $strand = shift;

	unless ($strand eq '+' || $strand eq '-') {
		confess "Error, strand value must be [+-]";
	}

	my $bit_val = ($strand eq '+') ? 0 : 1;
	$self->_set_bit_val(0x0020, $bit_val);
}

####
sub is_first_in_pair {
	my $self = shift;
	return($self->_get_bit_val(0x0040));
}

sub set_first_in_pair {
	my $self = shift;
	my $bit_val = shift;
	
	$self->_set_bit_val(0x0040, $bit_val);
	return;
}

####
sub is_second_in_pair {
	my $self = shift;
	return($self->_get_bit_val(0x0080));
}


sub set_second_in_pair {
	my $self = shift;
	my $bit_val = shift;

	$self->_set_bit_val(0x0080, $bit_val);
	return;
}



####
sub _get_bit_val {
	my $self = shift;
	my ($bit_position) = @_;

	my $flag = $self->get_flag();
	return($flag & $bit_position);
}


####
sub _set_bit_val {
	my $self = shift;
	my ($bit_position, $bit_val) = @_;

	unless (defined $bit_position && defined $bit_val) {
		confess "Error, need bit position and value";
	}
	
	my $flag = $self->get_flag();

	if ($bit_val) {
		$flag |= $bit_position;
	}
	else {
		# erase bit
		$flag &= ~$bit_position;
	}
	
	$self->set_flag($flag);
}


1; #EOM

