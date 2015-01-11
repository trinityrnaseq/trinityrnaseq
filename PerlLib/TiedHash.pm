#!/usr/local/bin/perl

package TiedHash;
use strict;
use warnings;
use DB_File;
use Carp;

=example

	my $tied_hash = new TiedHash( { create => "$pfam_db.inx" } );


    my $acc = "";

    while (<$fh>) {
	   chomp;
       my ($token, $rest) = split (/\s+/, $_, 2);
       if ($token eq 'NAME') {
	      $acc = $rest;
       }
       elsif ($token =~ /^(NC|TC|DESC|ACC)$/) {
	   my $key = "$acc$;$token";
       $tied_hash->store_key_value($key, $rest);
       print STDERR "storing: $key, $rest\n";
   }


=cut


sub new {
    my $packagename = shift;
    
    my $prefs_href = shift;
    
    if ($prefs_href && ! ref $prefs_href) {
        confess "Error, need hash reference with opts in constructor.\n";
    }
    

    my $self = { 
        index_filename => undef,
        tied_index => {},
        tie_invoked => 0,
    };
    
    bless ($self, $packagename);
    

    if (ref $prefs_href eq "HASH") {
        if (my $index_file = $prefs_href->{"create"}) {
            $self->create_index_file($index_file);
        }
        elsif ($index_file = $prefs_href->{"use"}) {
            $self->use_index_file($index_file);
        }
    }
            
    
    return ($self);
}

####
sub tie_invoked {
    my $self = shift;
    return ($self->{tie_invoked});
}


####
sub DESTROY {
    my $self = shift;
    if ($self->{index_filename}) {
        # hash must have been tied
        # so, untie it
        untie (%{$self->{tied_index}});
    }
}


####
sub create_index_file {
    my $self = shift;
    return ($self->make_index_file(@_));
}



####
sub make_index_file {
    my $self = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "need filename as parameter";
    }

    if (-e $filename) {
        unlink $filename or confess "cannot remove existing index filename $filename";
    }
    
    $self->{index_filename} = $filename;
    
    tie (%{$self->{tied_index}}, 'DB_File', $filename, O_CREAT|O_RDWR, 0666, $DB_BTREE);

    $self->{tie_invoked} = 1;
    
    return;
}


####
sub use_index_file {
    my $self = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "need filename as parameter";
    }
    
    unless (-s $filename) {
        confess "Error, cannot locate file: $filename\n";
    }
    
    $self->{index_filename} = $filename;
    
    tie (%{$self->{tied_index}}, 'DB_File', $filename, O_RDONLY, 0, $DB_BTREE);

    $self->{tie_invoked} = 1;

    #my @keys = $self->get_keys();
    #unless (@keys) {
    #    confess "Error, tried using $filename db, but couldn't perform retrievals.\n";
    #}
    
    return;

}


####
sub store_key_value {
    my ($self, $identifier, $value)  = @_;
    
    #my $num_keys = scalar ($self->get_keys());
    
    unless ($self->tie_invoked()) {
        confess "Error, cannot store key/value pair since tied hash not created.\n";
    }


    my $found = 0;
    while (! $found) {
        $self->{tied_index}->{$identifier} = $value;
        
        my $val = $self->get_value($identifier);
        if (defined $val) {
            $found = 1;
        }
        else {
            warn "Berkeley DB had trouble storing ($identifier); trying again.\n";
        }
    }
    
    return;
    
}


####
sub get_value {
    my $self = shift;
    my $identifier = shift;

      
    unless ($self->tie_invoked()) {
        confess "Error, cannot retrieve value from untied hash\n";
    }

    my $value = $self->{tied_index}->{$identifier};
    
    return ($value);
}


## 
sub get_keys {
    my $self = shift;
    
    unless ($self->tie_invoked()) {
        confess "Error, cannot retrieve values from untied hash\n";
    }
    
    return (keys %{$self->{tied_index}});
}


1; #EOM
