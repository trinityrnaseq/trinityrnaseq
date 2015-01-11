package Ascii_genome_illustrator;

use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	
	my ($molecule_name, $illustration_length) = @_;

	my $self = { name => $molecule_name,
				 illustration_length => $illustration_length,
				 features => [],
			 };

	bless ($self, $packagename);

	return ($self);
}


####
sub add_features {  ## accepts list of features
	my $self = shift;
	my @features = @_;

	foreach my $feature (@features) {
		unless (ref $feature eq 'ARRAY') { confess "Error, improper params; should be a list of feature array refs"; }
		$self->add_feature(@$feature);
	}


	return;
}
	
####
sub add_feature {
	my $self = shift;
	my ($feature_name, $feature_end5, $feature_end3, $glyph) = @_;

	unless ($feature_name && $feature_end5 =~ /^\d+$/ && $feature_end3 =~ /^\d+$/ && defined($glyph) ) { 
		confess "Error, improper params";
	}
	
	unless (length($glyph) == 1 && $glyph !~ /\s/) { confess "Error, glyph must be a single non ws character.";}
	
	my $orient = ($feature_end5 < $feature_end3) ? '+' : '-';
	
	my ($lend, $rend) = sort {$a<=>$b} ($feature_end5, $feature_end3);

	#print "Feature: $feature_name, $lend => $rend ($orient)\n";
	
	my $feature_struct = { name => $feature_name,
						   lend => $lend,
						   rend => $rend,
						   orient => $orient,
						   glyph => $glyph,
					   };
	
	push (@{$self->{features}}, $feature_struct);
	
	return;
}


####
sub get_features {
	my $self = shift;
	return (@{$self->{features}});
}



####
sub illustrate {
	my $self = shift;
	
	my ($mol_lend, $mol_rend) = @_;


	unless ($mol_lend && $mol_rend) { confess "invalid params"; }
	
	my $illustration_length = $self->{illustration_length};
	
	my $text = sprintf("%+20s  ", "$mol_lend-$mol_rend") . "[" . ("=" x ($illustration_length-2)) . "]\t" . $self->{name} . "\n";
	
	foreach my $feature ($self->get_features()) {
		
		my ($name, $lend, $rend, $orient, $glyph) = ($feature->{name}, $feature->{lend}, $feature->{rend}, $feature->{orient}, $feature->{glyph});
				
		my @illustration_array;
		## init to clean palette
		for (my $i = 0; $i < $illustration_length; $i++) { $illustration_array[$i] = " "; }
		
		my ($pos_lend, $pos_rend) = $self->_compute_palette_position([$lend, $rend], [$mol_lend, $mol_rend]);
		
		# draw
		for (my $i = $pos_lend; $i <= $pos_rend; $i++) { $illustration_array[$i] = $glyph; } 
		
		if ($orient eq '+') {
			$illustration_array[$pos_rend] = '>';
		}
		elsif ($orient eq '-')  {
			$illustration_array[$pos_lend] = '<';
		}
		else {
			confess "Don't recognize orient:$orient";
		}
		
	    $text .= sprintf ("%+20s$orient ", "$lend-$rend") . join ("", @illustration_array) . "\t$name\n";
	}
	
	return ($text);
}


####
sub _compute_palette_position {
	my $self = shift;
	
	my ($feature_coords_aref, $mol_coords_aref) = @_;

	my $illustration_length = $self->{illustration_length};
	
	my ($feature_lend, $feature_rend) = @$feature_coords_aref;
	my ($mol_lend, $mol_rend) = @$mol_coords_aref;

	unless ($feature_lend <= $mol_rend && $feature_rend >= $mol_lend) {
		## no overlap
		return (-1, -1);
	}

	if ($feature_lend < $mol_lend) {
		$feature_lend = $mol_lend;
	}

	if ($feature_rend > $mol_rend) {
		$feature_rend = $mol_rend;
	}

	my $mol_region_length = $mol_rend - $mol_lend + 1;

	my $delta_left = $feature_lend - $mol_lend + 1;
	my $delta_right = $feature_rend - $mol_lend + 1;

	my $pos_left = int($delta_left / $mol_region_length * $illustration_length + 0.5) - 1;
	$pos_left = 0 if $pos_left < 0;
	my $pos_right = int($delta_right / $mol_region_length * $illustration_length + 0.5) - 1;
	$pos_right = 0 if $pos_right < 0;
	
	return ($pos_left, $pos_right);
}



1; #EOM
