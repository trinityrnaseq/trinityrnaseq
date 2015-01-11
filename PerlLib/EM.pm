package EM;

use strict;
use warnings;
use SAM_entry;
use Carp;

use vars qw($TOKEN);
BEGIN {
	$TOKEN = "$;";
}

my $MIN_EFF_TRANS_LENGTH = 10; # // large enough to keep abundance estimates from going absurdly high for short transcripts.

####
sub  new {
	my $packagename = shift;
	my ($transcript_seqs_href, $fragment_length) = @_;
	
    unless ($fragment_length) {
        $fragment_length = 1; # negligible with respect to the length of the transcript
    }


	# transcript_seqs_href:  ( seq_acc => sequence, ...)
	# trans_reads_aref: ( [trans_accA, read_accA], [trans_accB, read_accB], ...)
	
	unless ($transcript_seqs_href) {
		confess "Error, need param: transcript_seqs_href";
	}
	
	

	my $self = {
		
		# mapping info
		_trans_to_multi_map_counts => {}, # trans => { transA\$;transB\$;transC... } => count_of_reads

		# model params
		_ENt => {}, # trans => expected read count
		_theta => {}, # trans => fraction of all reads
		
		# transcript info
		_trans_lengths => {}, # trans => length
        _frag_length => $fragment_length,
        
				
		# misc
		_total_reads => 0, # total read count

	};
	

	bless ($self, $packagename);

	$self->_init_trans_lengths($transcript_seqs_href);
	
	return($self);

}


####
sub _init_trans_lengths {
	my $self = shift;
	my ($transcript_seqs_href) = @_;
	
	foreach my $transcript (keys %$transcript_seqs_href) {
		
		my $sequence = $transcript_seqs_href->{$transcript} or confess "Error, no sequence for transcript: $transcript";
		
		$self->{_trans_lengths}->{$transcript} = length($sequence);
		
	}

	return;
}

sub run {
	my $self = shift;
	my (%settings) = @_;

	my $max_iterations = $settings{max_iterations} || 1000;
	my $min_delta_ML = $settings{min_delta_ML} || 0.01;
    my $verbose_flag = $settings{verbose} || 0;

	
	$self->_init_theta();
	
	#########################
	## Do the EM:
	
	my $prev_ML_val = undef;
	
	my $round = 0;
	my $delta = 1;
	
	print STDERR "\tEM started\n";
	while (1) {
	
		$round++;
		
		$self->_E_step();

		$self->_M_step();

		
		my $ml = $self->_compute_Max_Likelihood();
		
		if (defined $prev_ML_val) {
			$delta = $ml - $prev_ML_val;
		}
		
		print STDERR "\rML[$round]: $ml       $delta       ";
		
		$prev_ML_val = $ml;
		
        if ($verbose_flag) {
            print "\n\nEM_results, round: $round:\n";
            $self->report_results();
            print "\n";
        }
        
		if ($round >= $max_iterations || $delta  < $min_delta_ML) {
			
			last;
		}

	}
	print STDERR "\n\tEM finished.\n";
	
	return;
}


####
sub _init_theta {
	my $self = shift;

	###########################
	## init theta:
	##   assume mutliply mapped reads equally divided among their mapped locations for init
	
	my @transcripts = $self->_get_all_transcripts();

	my $total_reads = $self->{_total_reads};
	
	foreach my $transcript (@transcripts) {
		
		my $read_count_sum = 0;
		
		my @combos = $self->_get_multi_map_list_including_transcript($transcript);
		
		foreach my $combo (@combos) {
			
			my @other_trans = $self->_get_transcripts_in_combo($combo);
			
			my $num_other = scalar(@other_trans);
			my $count = $self->_get_multi_map_read_count_for_combo($transcript, $combo);
			
			$read_count_sum += ($count/$num_other);
		}
				
		$self->{_theta}->{$transcript} = $read_count_sum / $total_reads;

	}

	return;
}

####
sub _get_all_transcripts {
	my $self = shift;
	
	my @transcripts = keys %{$self->{_trans_to_multi_map_counts}};
	
	return(@transcripts);
}


####
sub _get_reads_mapped_to_transcript {
	my $self = shift;
	my ($transcript) = @_;

	if (exists $self->{_trans_to_reads_href}->{$transcript}) {
		my @reads = @{$self->{_trans_to_reads_href}->{$transcript}};
		return(@reads);
	}
	else {
		confess "Error, no reads mapped to transcript $transcript";
	}
}

####
sub _get_transcript_length {
	my $self = shift;
	my ($transcript) = @_;

	return($self->{_trans_lengths}->{$transcript});
}




####
sub _E_step {
	my $self = shift;

	my @transcripts = $self->_get_all_transcripts();
	
	## do the E
	foreach my $transcript (@transcripts) {
		
		my $expected_read_count = 0;
		
		
		
		my $trans_length = $self->_get_transcript_length($transcript) or die "Error, no length for transcript: $transcript";
        my $eff_trans_length = $trans_length - $self->{_frag_length} + 1;
        if ($eff_trans_length < $MIN_EFF_TRANS_LENGTH) {
            $eff_trans_length = $MIN_EFF_TRANS_LENGTH;
        }
        

		my $theta_t = $self->{_theta}->{$transcript};
		if (! defined $theta_t) {
			die "Error, no theta($transcript), $theta_t";
		}
		
		
		my @combos = $self->_get_multi_map_list_including_transcript($transcript);
		
		#print "Combos: @combos\n";
		

	    foreach my $combo (@combos) {
						
			my $combo_count = $self->_get_multi_map_read_count_for_combo($transcript, $combo);
			
			#print "Got combo: $combo, count: $combo_count\n";
						
			my @other_transcripts = split(/$TOKEN/, $combo);
			if (scalar @other_transcripts == 1) {
				# only current transcript:
				unless ($other_transcripts[0] eq $transcript) {
					confess "Error, mapped single transcript should be the current transcript";
				}
				$expected_read_count += $combo_count;
			}
			else {
				## compute partial mapping.
				
							
				my $numerator = (1 / $eff_trans_length) * $theta_t;
                #my $numerator = $theta_t;
                
			
				# demonator: sum for all t: P(r|t) * theta(t)
				my $denominator = 0;
			
				
				#print STDERR "Other transcripts: @other_transcripts\n";
				
				foreach my $other_transcript (@other_transcripts) {  # other transcripts includes the current transcript as well.
				
					my $other_trans_len = $self->_get_transcript_length($other_transcript);
					unless ($other_trans_len) {
						die "Error, no trans length for $other_transcript";
					}
					
                    my $eff_other_trans_length = $other_trans_len - $self->{_frag_length} + 1;
                    if ($eff_other_trans_length < $MIN_EFF_TRANS_LENGTH) {
                        $eff_other_trans_length = $MIN_EFF_TRANS_LENGTH;
                    }
                    
					my $other_theta = $self->{_theta}->{$other_transcript};
					
					if (! defined $other_theta) {
						die "no theta for $other_transcript, $other_theta";
					}
					
					my $val = (1 / $eff_other_trans_length) * $other_theta;
                    #my $val = $other_theta;
                    
					$denominator += $val;
				}
			
				#print "Fractional read count for $transcript = $numerator / $denominator\n";
			
				my $fractional_read_count = $combo_count * ($numerator/$denominator);
				
				
				$expected_read_count += $fractional_read_count;
			}
		}
		
		$self->{_ENt}->{$transcript} = $expected_read_count;
	}
	
	return;
}
						
			
####
sub _M_step {
	my $self = shift;
	
	my $sum_NT = 0;

	foreach my $E (values %{$self->{_ENt}}) {
		$sum_NT += $E;
	}
	
	## theta = ENt(t) / sum( for all t: ENt(t) )

	
	foreach my $transcript ($self->_get_all_transcripts()) {
		my $new_theta = $self->{_ENt}->{$transcript} / $sum_NT;
		
		$self->{_theta}->{$transcript} = $new_theta;
	}

	return;
}


####
sub _compute_Max_Likelihood {
	my $self = shift;
	
	# ML = sum:t (ENt * log(theta(t)))

	my $ML = 0;

	foreach my $transcript ($self->_get_all_transcripts()) {
		
		my $theta = $self->{_theta}->{$transcript};
		
		if ($theta > 0) {
			
			my $ENt = $self->{_ENt}->{$transcript};

			$ML += $ENt * log($theta);
			
		}
	}

	return($ML);
}


####
sub get_results {
	my $self = shift;
	
	my $total_reads = $self->{_total_reads};

    my @results;
	
	foreach my $transcript (sort $self->_get_all_transcripts()) {
		
		my $num_frags = $self->{_ENt}->{$transcript};
		my $len = $self->_get_transcript_length($transcript);
		
		$num_frags = sprintf("%.2f", $num_frags);
		
		my $fpkm = $num_frags / ($len/1e3) / ($total_reads/1e6);
		$fpkm = sprintf("%.2f", $fpkm);

		my ($num_unique_reads, $num_multi_map_reads) = $self->count_reads_mapped_to_transcript($transcript);
		
        my $struct = { trans_id => $transcript,
                       length => $len,
                       unique_map => $num_unique_reads,
                       multi_map => $num_multi_map_reads,
                       expected_map => $num_frags,
                       FPKM => $fpkm,
        };

        push (@results, $struct);


		
	}
	

	return (@results);
}


sub print_results {
    my $self = shift;
    
    print $self->report_results();
    return;
}



sub report_results {
    my ($self) = shift;

    my @results = $self->get_results();
        
    my $out_text = join("\t", "#transcript", "trans_length", "unique_map", "multi_map", "EM_frag_count", "FPKM") . "\n";
    
    foreach my $result (@results) {
            
        $out_text .= join("\t", 
                          $result->{trans_id},
                          $result->{length},
                          $result->{unique_map},
                          $result->{multi_map},
                          $result->{expected_map},
                          $result->{FPKM}) . "\n";
    }
    
    $out_text .= "\n"; # last spacer
    
    return ($out_text);

}






####
sub _get_multi_map_list_including_transcript {
	my $self = shift;
	my ($transcript) = @_;

	my @combos = keys %{$self->{_trans_to_multi_map_counts}->{$transcript}};
	return(@combos);
}

####
sub _get_multi_map_read_count_for_combo {
	my $self = shift;
	my ($transcript, $combo) = @_;
	
	my $val = $self->{_trans_to_multi_map_counts}->{$transcript}->{$combo};
	
	return($val);
}

####
sub _get_transcripts_in_combo {
	my $self = shift;
	my ($combo) = @_;

	my @transcripts = split(/$TOKEN/, $combo);
	return(@transcripts);
}


####
sub _create_combo {
	my $self = shift;
	my @transcripts = @_;

	my $combo = join($TOKEN, sort @transcripts);
	
	return($combo);
}

####
sub add_read { # really add_fragment: count pairs only once.
	my $self = shift;
	my @transcripts = @_;

	my $combo = $self->_create_combo(@transcripts);

	foreach my $transcript (@transcripts) {
		$self->{_trans_to_multi_map_counts}->{$transcript}->{$combo}++;
	}

	$self->{_total_reads}++;

	return;
}

####
sub count_reads_mapped_to_transcript {
	my $self = shift;
	my ($transcript) = @_;
	
	my $unique_read_count = 0;
	my $multi_map_count = 0;


	my @combos = $self->_get_multi_map_list_including_transcript($transcript);                                                       
	
	foreach my $combo (@combos) {                                                                                                    
                                                                                                                                         
		my $combo_count = $self->_get_multi_map_read_count_for_combo($transcript, $combo);                                           
		
		
		my @other_transcripts = $self->_get_transcripts_in_combo($combo);
		
		if (scalar @other_transcripts == 1) {                                                                                        
			# only current transcript:                                                                                               
			$unique_read_count += $combo_count;
		}                                                                                                                            
		else {                                   
			$multi_map_count += $combo_count;
		}
	}

	return($unique_read_count, $multi_map_count);
}

1; #EOM

