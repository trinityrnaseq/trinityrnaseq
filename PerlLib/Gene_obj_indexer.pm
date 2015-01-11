#!/usr/local/bin/perl

package Gene_obj_indexer;
use strict;
use warnings;
use base qw(TiedHash);
use Gene_obj;
use Storable qw (thaw nfreeze);
use Carp;


####
sub new {
    my $packagename = shift;
    
    my $self = $packagename->SUPER::new(@_);
 
    return ($self);
   
}

####
sub store_gene {
    my ($self, $identifier, $gene_obj)  = @_;
    

    unless (ref $gene_obj) {
        confess "Error, no gene_obj as param";
    }
    
    my $blob = nfreeze ($gene_obj);
    
    my $success = 0;
    
    while (! $success) {
        $self->store_key_value($identifier, $blob);
        
        eval {
            my $gene_obj = $self->get_gene($identifier);
            
        };
        if ($@) {
            warn "error trying to store gene $identifier using berkeley db.  Trying again...\n";
        }
        else {
            # worked.
            $success = 1;
        }
    }
    
}


####
sub get_gene {
    my $self = shift;
    my $identifier = shift;

    my $blob = $self->get_value($identifier);
    
    unless ($blob) {
        confess "Error, no gene obj retrieved based on identifier $identifier";
    }

    my $gene_obj = thaw($blob);
    unless (ref $gene_obj) {
        confess "Error retrieving gene_obj based on identifier $identifier.  Data retrieved but not thawed properly.\n";
    }
    
    return ($gene_obj);
}


1; #EOM

