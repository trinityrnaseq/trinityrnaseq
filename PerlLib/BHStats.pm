package BHStats;

require Exporter;
our @ISA = qw (Exporter);
our @EXPORT = qw (binomial_probability_sum_k_to_n
				  binomial_probability_sum_k_to_0
				  binomial_probability
				  binomial_coefficient
				  factorial
				  stDev
				  stdErr
				  median
				  avg
				  CorrelationCoeff
				  geometric_mean
				  min
				  max
				  sum
                  tukey_biweight
				  );


use strict;

sub binomial_probability_sum_k_to_n {
    my ($n,$k,$p) = @_;
    my $sum = 0;
    for (my $i = $k; $i <= $n; $i++) {
		my $binProb = binomial_probability($n,$i,$p);
		$sum += $binProb;
    }
    return ($sum);
}

sub binomial_probability_sum_k_to_0 {
    my ($n,$k,$p) = @_;
    my $sum = 0;
    for (my $i = $k; $i >= 0; $i--) {
		my $binProb = binomial_probability($n,$i,$p);
		$sum += $binProb;
    }
    return ($sum);
}




sub binomial_probability {
    my ($n_observations, $k_successes, $p_probability) = @_;
    
    my ($n, $k, $p) = ($n_observations, $k_successes, $p_probability);
    
    ### Given B(n,p), find P(X=k)
    
    my $binomial_prob = binomial_coefficient($n,$k) * ($p**$k) * (1-$p)**($n-$k);
	
    return ($binomial_prob);
    
}


sub binomial_coefficient {
    my ($n_things, $k_at_a_time) = @_;
    
    my $number_of_k_arrangements = (factorial($n_things)) / ( factorial($k_at_a_time) * factorial($n_things-$k_at_a_time) );

    return ($number_of_k_arrangements);
}


sub factorial {
    my $x = shift;
    $x = int($x);
    my $factorial = 1;
    while ($x > 1) {
		$factorial *= $x;
		$x--;
    }
    return ($factorial);
}


sub stDev {
    # standard deviation calculation
    my @nums = @_;
    @nums = sort {$a<=>$b} @nums;
	
	
    my $avg = avg(@nums);
    my $count_eles = scalar(@nums);
    
    ## sum up the sqr of diff from avg
    my $sum_avg_diffs_sqr = 0;
    foreach my $num (@nums) {
        my $diff = $num - $avg;
        my $sqr = $diff**2;
        $sum_avg_diffs_sqr += $sqr;
    }
    my $stdev = sqrt ($sum_avg_diffs_sqr/($count_eles-1));
    return ($stdev);
}

####
sub stdErr {
	my @vals = @_;

	my $stdev = &stDev(@vals);
	
	my $num_vals = scalar(@vals);

	my $stdErr = $stdev / sqrt($num_vals);

	return($stdErr);
}


sub median {
    my @nums = @_;
    
	@nums = sort {$a<=>$b} @nums;
	
    my $count = scalar (@nums);
    if ($count %2 == 0) {
        ## even number:
        my $half = $count / 2;
        return (avg ($nums[$half-1], $nums[$half]));
    }
    else {
        ## odd number. Return middle value
        my $middle_index = int($count/2);
        return ($nums[$middle_index]);
    }
}

sub avg {
    my @nums = @_;
    my $total = $#nums + 1;
    my $sum = 0;
    foreach my $num (@nums) {
		$sum += $num;
    }
    my $avg = $sum/$total;
    return ($avg);
}


sub CorrelationCoeff {
    my ($x_aref, $y_aref) = @_;
    my @x = @$x_aref;
    my @y = @$y_aref;
    
    my $total = $#x + 1;
    my $avg_x = avg(@x);
    my $avg_y = avg(@y);
    
    my $stdev_x = stDev(@x);
    my $stdev_y = stDev(@y);
    
    # sum part of equation
    my $summation = 0;
    for (my $i = 0; $i < $total; $i++) {
		my $x_val = $x[$i];
		my $y_val = $y[$i];
		
		my $x_part = ($x_val - $avg_x)/$stdev_x;
		my $y_part = ($y_val - $avg_y)/$stdev_y;
		
		$summation += ($x_part * $y_part);
    }
    
    my $cor = (1/($total-1)) * $summation;
    
    return ($cor);
}


####
sub geometric_mean {
    my @entries = @_;
    
    my $num_entries = scalar (@entries);
    unless ($num_entries) {
		return (undef);
    }
    
    ## All entries must be > 0
    my $logsum = 0;
    foreach my $entry (@entries) {
		unless ($entry > 0) {
			return (undef);
		}
		$logsum += log ($entry);
    }
    
    my $geo_mean = exp ( (1/$num_entries) * $logsum);
    
    return ($geo_mean);
}


####
sub min {
	my @vals = @_;
	
	@vals = sort {$a<=>$b} @vals;
	
	my $min_val = shift @vals;
	
	return ($min_val);
}

####
sub max {
	my @vals = @_;
	
	@vals = sort {$a<=>$b} @vals;
	
	my $max_val = pop @vals;
	
	return ($max_val);
}


####
sub sum {
	my @vals = @_;

	my $x = 0;
	foreach my $val (@vals) {
		$x += $val;
	}
	
	return ($x);
}


=Rcode

> tukey.biweight
function (x, c = 5, epsilon = 1e-04) 
{
    m <- median(x)
    s <- median(abs(x - m))
    u <- (x - m)/(c * s + epsilon)
    w <- rep(0, length(x))
    i <- abs(u) <= 1
    w[i] <- ((1 - u^2)^2)[i]
    t.bi <- sum(w * x)/sum(w)
    return(t.bi)
}
=cut

####
sub tukey_biweight {
    my (@x) = @_;

    my $m = median(@x);
    
    my $s;
    {
        my @y;
        foreach my $val (@x) {
            my $t = abs($val - $m);
            push (@y, $t);
        }
        $s = median(@y);
    }

        
    my @u;
    {
        my $epsilon = 1e-4;
        my $c = 5;
                
        foreach my $val (@x) {
            my $t = $val - $m;
            $t /= ($c * $s + $epsilon);
            push (@u, $t);
        }
    }
    
    my @i;
    {
        foreach my $val (@u) {
            my $t = (abs($val) <= 1) ? 1:0;
            push (@i, $t);
        }
    }

    my @w;
    {
        foreach my $val (@u) {
            my $i = shift @i;
            my $t = ( (1 - $val**2) **2) * $i;
            push (@w, $t);
        }
    }
                
    my $bi;
    {
        my $sum_w = 0;
        for (my $i = 0; $i < @w; $i++) {
            my $w = $w[$i];
            my $x = $x[$i];
            $bi += $w * $x;
            $sum_w += $w;
        }
        $bi /= $sum_w;
    }

    return($bi);
}




1; #EOM
