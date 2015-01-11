package Ktree;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    
    my $self = { _root => KtreeNode->new("", 0) };

    bless ($self, $packagename);

    return($self);
}

sub add_kmer {
    my $self = shift;
    my ($kmer) = @_;


    unless (defined $kmer) {
        confess "error, require param kmer";
    }

    my $root_node = $self->{_root};

    my @seq = split(//, $kmer);

    my $node = $root_node; 
    do {
        my $char = shift @seq;
        $node = $node->get_child($char);
    } while (@seq);
    
    $node->set_val( $node->get_val() + 1 );
    
    return;
}

sub report_kmer_counts {
    my $self = shift;
    
    my $root_node = $self->{_root};
    
    &_recurse_through_kmer_counts("", $root_node);
    
    return;
}

sub _recurse_through_kmer_counts {
    my ($prefix, $node) = @_;

    my $char = $node->get_char();
    
    my @children_chars = $node->get_children_chars();
    
    if (@children_chars) {
        foreach my $child_char (@children_chars) {
            my $child_node = $node->get_child($child_char);
            &_recurse_through_kmer_counts($prefix . $char, $child_node);
        }
    }
    else {
        # base case
        my $val = $node->get_val();
        print join("\t", $prefix . $char, $val) . "\n";
    }

    return;
}



package KtreeNode;

use strict;
use warnings;
use Carp;


sub new {
    my $packagename = shift;
    my ($char, $val) = @_;

    unless (defined $char && defined $val) {
        confess "Error, require (character, val) as parameter";
    }

    my $self = { char => $char,
                 val => $val,
                 children => {},
    };
    
    bless ($self, $packagename);
    
    return($self);
}




####
sub get_child {
    my $self = shift;
    my ($char) = @_;

    unless (defined $char) {
        confess "error, parameter 'char' required";
    }

    my $child = $self->{children}->{$char};
    unless (ref $child) {
        
        $child = $self->{children}->{$char} = new KtreeNode($char, 0);
    }

    return($child);
}


sub get_children_chars {
    my $self = shift;
    
    my @chars = keys %{$self->{children}};
    return(@chars);
}


sub get_char {
    my $self = shift;
    return($self->{char});
}


####
sub get_val {
    my $self = shift;
    return($self->{val});
}

####
sub set_val {
    my $self = shift;
    my $val = shift;
    unless (defined $val) {
        confess "error, require val as param";
    }

    $self->{val} = $val;
    
    return;
}


1; #EOM
