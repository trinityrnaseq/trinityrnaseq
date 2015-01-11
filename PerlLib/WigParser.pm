package WigParser;

use strict;
use warnings;
use Carp;

####
sub new {
    my $packagename = shift;
    my ($wig_file) = @_;

    unless ($wig_file) {
        confess "Error, need wig_filename as parameter to constructor";
    }

    my $self = { wig_file => $wig_file,
                 wig_file_fh => undef,
                 contig_to_seek_pos_href => {},
             };

    bless ($self, $packagename);
    
    $self->_init_wig_info();


    return($self);
}

####
sub _init_wig_info {
    my $self = shift;

    my $wig_file = $self->{wig_file};
    

    open (my $fh, $wig_file) or die "Error, cannot open file $wig_file";
    $self->{wig_file_fh} = $fh;
    
    while (<$fh>) {
        if (/chrom=(\S+)/) {
			my $scaff = $1;
            my $filepos = tell($fh);
            $self->{contig_to_seek_pos_href}->{$scaff} = $filepos;
        }
    }
    

    return;
}

####
sub get_contig_list {
    my $self = shift;
    
    my @contigs = keys %{$self->{contig_to_seek_pos_href}};

    return(@contigs);
}

####
sub get_wig_array {
    my $self = shift;
    my ($contig, $extended_flag) = @_;
    
    
    my $seekpos = $self->{contig_to_seek_pos_href}->{$contig};
    
    unless (defined $seekpos) {
        confess "Error, cannot find seek position entry for contig: $contig";
    }

    my $fh = $self->{wig_file_fh};

    seek($fh, $seekpos, 0); 

    #print STDERR "-retrieving wig array for $contig at pos: $seekpos\n";
    
    my @wig_array;

    while (<$fh>) {
        
        if (/chrom=/) { last; }
        
        chomp;
        my ($pos, $val, @rest) = split(/\s+/);
        if ($pos > 1e12) {
            confess "Error, position $pos is out of range of max contig position value set at 1e12";
        }
        if ($extended_flag) {
            $wig_array[$pos] = [$val, @rest];
        }
        else {
            $wig_array[$pos] = $val;
        }
        
    }
    
    # fill in missing entries
    foreach my $val (@wig_array) {
        if (! defined ($val)) {
            $val = 0;
        }
    }


    return(@wig_array);
}



#####################################################################
## Static method that uses LOTS of memory (original implementation)

sub parse_wig {
	my ($wig_file) = @_;

	my %scaff_to_coverage;
    
	print STDERR "-retrieving max positions per scaffold\n";
	my %max_vals_for_scaffs = &_parse_max_scaff_lengths($wig_file);
	
	print STDERR "-preallocating memory for coverage\n";
	## preallocate arrays for coverage

	my $sum_pos = 0;
	foreach my $scaff (keys %max_vals_for_scaffs) {
		my $max_val = $max_vals_for_scaffs{$scaff};

		my @cov_vals = ();
		$#cov_vals = $max_val;
		
		for my $index (0..$max_val) {
			$cov_vals[$index] = int(0);
		}
		$scaff_to_coverage{$scaff} = \@cov_vals;
		
		$sum_pos += $max_val;
	
	}

	print STDERR "- $sum_pos bases represented\n";

	print STDERR "-populating coverage data into memory.\n";
	
	my $scaff = undef;

	open (my $fh, $wig_file) or die "Error, cannot open file $wig_file";
	while (<$fh>) {
		if (/^track/) { next; };
		if (/chrom=(\S+)/) {
			$scaff = $1;
			next;
		}
		chomp;
		my ($pos, $val) = split(/\s+/);
		$scaff_to_coverage{$scaff}->[$pos] = $val;
		
	}
	
	close $fh;

	return(%scaff_to_coverage);
}



### Private

sub _parse_max_scaff_lengths {
	my ($wig_file) = @_;

	my %scaff_to_max_vals;
	
	my $scaff = "";


	my $counter = 0;
	open (my $fh, $wig_file) or die "Error, cannot open file $wig_file";
	while (<$fh>) {
		if (/^track/) { next; };
		if (/chrom=(\S+)/) {
			$scaff = $1;
			next;
		}
		chomp;
		my ($pos, $val) = split(/\s+/);
		$scaff_to_max_vals{$scaff} = $pos;
		
		#print "scaff: $scaff\t$pos=> $val\n";
		
		$counter++;

		#if ($counter > 100000) { last; }
	}
	
	close $fh;

	return(%scaff_to_max_vals);
}

1; #EOM

