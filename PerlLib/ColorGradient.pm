package ColorGradient;

use strict;
use warnings;


# RGB values off (0), on (1), increasing (I), or decreasing (D).
my @ColorPhases = ( ['1', 'I', '0'],
					['D', '1', '0'],
					['0', '1', 'I'],
					['0', 'D', '1'],
					['I', '0', '1'],
					);

my $num_color_phases = scalar @ColorPhases;
my $discrete_color_phase_percent = 1 /$num_color_phases * 100;

####
sub get_RGB_gradient {
  my ($num_colors) = @_;
  
  ## Returns a list of [R,G,B], ... values.

  unless ($num_colors =~ /^\d+$/) {
	die "Error, $num_colors should be a number";
  }

  my @colors;
  $num_colors--;
  for my $color_entry (0..$num_colors) {
	
	my $percentage = ($color_entry-0.000001) / $num_colors * 100; # don't ever want to reach 100%
	
	my $index = int ($percentage / $discrete_color_phase_percent);

	my $color_phase = $ColorPhases[$index];

	my $ratio_into_phase = ($percentage - $discrete_color_phase_percent*$index) / $discrete_color_phase_percent;
	
	my $rgb_color = &_get_color($color_phase, $ratio_into_phase);

	push (@colors, $rgb_color);
  }
  
  return (@colors);
}

###
sub convert_RGB_hex {
  my @rgb_aref_vals = @_;
  
  my @hex_vals;
  foreach my $rgb_aref (@rgb_aref_vals) {
	my ($r, $g, $b) = @$rgb_aref;
	
	my $hex_r = sprintf("%02x", $r);
	my $hex_g = sprintf("%02x", $g);
	my $hex_b = sprintf ("%02x", $b);
	
	push (@hex_vals, '#' . $hex_r . $hex_g . $hex_b);
  }
  return (@hex_vals);
}

####
sub _get_color {
  my ($color_phase, $ratio_into_phase) = @_;

  my $variable_rgb_val = int($ratio_into_phase * 255 + 4/9);
  
  my @phase_vals = @$color_phase;
  

  my @ret_rgb;
  foreach my $phase_val (@phase_vals) {
	my $rgb_val;
	if ($phase_val eq '1') {
	  $rgb_val = 255;
	}
	elsif ($phase_val eq '0') {
	  $rgb_val = 0;
	}
	elsif ($phase_val eq 'I') {
	  $rgb_val = $variable_rgb_val;
	}
	elsif ($phase_val eq 'D') {
	  $rgb_val = 255 - $variable_rgb_val;
	}
	else {
	  die "Error, unrecognized phase value of $phase_val";
	}
	
	push (@ret_rgb, $rgb_val);
  }

  return ([@ret_rgb]);
}




1;
	

