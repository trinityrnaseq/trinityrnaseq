package main;
our $SEE = 0;


package CDNA::Overlap_assembler;

use strict;

sub new {
    my $packagename = shift;
    my $self = {
	node_list => []
	};
    bless ($self, $packagename);
    return ($self);
}

sub add_cDNA {
    my $self = shift;
    my ($cdna_acc, $end5, $end3) = @_;
    my ($cdna_start, $cdna_stop) = sort {$a<=>$b} ($end5, $end3);
    my $node = Cdna_node->new($cdna_acc, $cdna_start, $cdna_stop);
    push (@{$self->{node_list}}, $node);
}


####
sub build_clusters {
    my $self = shift;
    my $node_list_aref = $self->{node_list};
    @{$node_list_aref} = sort {$a->{lend}<=>$b->{lend}} @{$node_list_aref}; #sort by lend coord.
    ## set indices
    for (my $i = 0; $i <= $#{$node_list_aref}; $i++) {
	$node_list_aref->[$i]->{myIndex} = $i;
    }
    
    my @clusters;
    my $first_node = $node_list_aref->[0];
    my $start_pos = 0;
    my ($exp_left, $exp_right) = ($first_node->{lend}, $first_node->{rend});
    print $first_node->{acc} . " ($exp_left, $exp_right)\n" if $SEE;
    for (my $i = 1; $i <= $#{$node_list_aref}; $i++) {
	my $curr_node = $node_list_aref->[$i];
	my ($lend, $rend) = ($curr_node->{lend}, $curr_node->{rend});
	print $curr_node->{acc} . " ($lend, $rend)\n" if $SEE;
	if ($exp_left <= $rend && $exp_right >= $lend) { #overlap
	    $exp_left = &min($exp_left, $lend);
	    $exp_right = &max($exp_right, $rend);
	    print "overlap. New expanded coords: ($exp_left, $exp_right)\n" if $SEE;
	} else {
	    print "No overlap; Creating cluster: " if $SEE;
	    my @cluster;
	    for (my $j=$start_pos; $j < $i; $j++) {
		my $acc = $node_list_aref->[$j]->{acc};
		push (@cluster, $acc);
		print "$acc, " if $SEE;
	    }
	    push (@clusters, [@cluster]);
	    $start_pos = $i;
	    ($exp_left, $exp_right) = ($lend, $rend);
	    print "\nResetting expanded coords: ($lend, $rend)\n" if $SEE;
	}
    }

    print "# Adding final cluster.\n" if $SEE;
    if ($start_pos != $#{$node_list_aref}) {
	print "final cluster: " if $SEE;
	my @cluster;
	for (my $j = $start_pos; $j <= $#{$node_list_aref}; $j++) {
	   my $acc = $node_list_aref->[$j]->{acc};
	   print "$acc, " if $SEE;
	   push (@cluster, $acc);
       }
	push (@clusters, [@cluster]);
	print "\n" if $SEE;
    } else {
	my $acc = $node_list_aref->[$start_pos]->{acc};
	push (@clusters, [$acc]);
	print "adding final $acc.\n" if $SEE;
    }
    return (@clusters);
}

sub min {
    my (@x) = @_;
    @x = sort {$a<=>$b} @x;
    my $min = shift @x;
    return ($min);
}

sub max {
    my @x = @_;
    @x = sort {$a<=>$b} @x;
    my $max = pop @x;
    return ($max);
}

#################################################################
package Cdna_node;
use strict;

sub new {
    my $packagename = shift;
    my ($acc, $lend, $rend) = @_;
    my $self = { acc=>$acc,
		 lend=>$lend,
		 rend=>$rend,
		 myIndex=>undef(),
		 overlapping_indices=>[]
		 };
    bless ($self, $packagename);
    return ($self);
}

1; #EOM




