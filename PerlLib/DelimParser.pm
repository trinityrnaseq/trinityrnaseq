#!/usr/bin/env perl

# classes for DelimParser::Reader and DelimParser::Writer

package DelimParser;
use strict;
use warnings;
use Carp;

####
sub new {
    my ($packagename, $fh, $delimiter) = @_;
    
    unless ($fh && $delimiter) {
        confess "Error, need filehandle and delimiter params";
    }
    

    my $self = {  _delim => $delimiter,
                  _fh => $fh,

                  # set below in _init()
                  _column_headers => [],
    };

    
    bless ($self, $packagename);
        
    return($self);
}


####
sub get_fh {
    my $self = shift;
    return($self->{_fh});
}

####
sub get_delim {
    my $self = shift;
    return($self->{_delim});
}

####
sub get_column_headers {
    my $self = shift;
    return(@{$self->{_column_headers}});
}

####
sub set_column_headers {
    my $self = shift;
    my (@columns) = @_;

    $self->{_column_headers} = \@columns;
    
    return;
}

####
sub get_num_columns {
    my $self = shift;
    return(length($self->get_column_headers()));
}


###
sub reconstruct_header_line {
    my $self = shift;
    my @column_headers = $self->get_column_headers();

    my $header_line = join("\t", @column_headers);
    return($header_line);
}

###
sub reconstruct_line_from_row {
    my $self = shift;
    my $row_href = shift;
    unless ($row_href && ref $row_href) {
        confess "Error, must set row_href as param";
    }

    my @column_headers = $self->get_column_headers();

    my @vals;
    foreach my $col_header (@column_headers) {
        my $val = $row_href->{$col_header};
        push (@vals, $val);
    }

    my $row_text = join("\t", @vals);

    return($row_text);
        
}


##################################################
package DelimParser::Reader;
use strict;
use warnings;
use Carp;
use Data::Dumper;

our @ISA;
push (@ISA, 'DelimParser');

sub new {
    my ($packagename) = shift;
    my $self = $packagename->DelimParser::new(@_);
    
    $self->_init();
    
    return($self);
}


####
sub _init {
    my $self = shift;
    
    my $fh = $self->get_fh();
    my $delim = $self->get_delim();

    my $header_row = <$fh>;
    chomp $header_row;

    unless ($header_row) {
        confess "Error, no header row read.";
    }
    
    my @fields = split(/$delim/, $header_row);
        
    $self->set_column_headers(@fields);
    
    
    return;
}

####
sub get_row {
    my $self = shift;
    
    my $fh = $self->get_fh();
    my $line = <$fh>;
    unless ($line) {
        return(undef); # eof
    }
    
    chomp $line;
    
    my $delim = $self->get_delim();
    my @fields = split(/$delim/, $line);
    
    my @column_headers = $self->get_column_headers();

    my $num_col = scalar (@column_headers);
    my $num_fields = scalar(@fields);

    if ($num_col != $num_fields) {
        confess "Error, line is lacking $num_col fields: " . Dumper(\@column_headers);
    }
    
    my %dict;
    foreach my $colname (@column_headers) {
        my $field = shift @fields;
        $dict{$colname} = $field;
    }
    
    return(\%dict);
}


##################################################

package DelimParser::Writer;
use strict;
use warnings;
use Carp;

our @ISA;
push (@ISA, 'DelimParser');

sub new {
    my ($packagename) = shift;
    my ($ofh, $delim, $column_fields_aref, $FLAGS) = @_;
        
    ## FLAGS can be:
    #  NO_WRITE_HEADER|...

    unless (ref $column_fields_aref eq 'ARRAY') {
        confess "Error, need constructor params: ofh, delim, column_fields_aref";
    }
    
    my $self = $packagename->DelimParser::new($ofh, $delim);
 
    $self->_initialize($column_fields_aref, $FLAGS);
    
    return($self);
}


####
sub _initialize {
    my $self = shift;
    my $column_fields_aref = shift;
    my $FLAGS = shift;
    
    unless (ref $column_fields_aref eq 'ARRAY') {
        confess "Error, require column_fields_aref as param";
    }
    
        
    my $ofh = $self->get_fh();
    my $delim = $self->get_delim();

    
    
    $self->set_column_headers(@$column_fields_aref);
    
    unless (defined($FLAGS) && $FLAGS =~ /NO_WRITE_HEADER/) {
        my $output_line = join($delim, @$column_fields_aref);
        print $ofh "$output_line\n";
    }
    
    
    return;
}


####
sub write_row {
    my $self = shift;
    my $dict_href = shift;

    unless (ref $dict_href eq "HASH") {
        confess "Error, need dict_href as param";
    }

    my $num_dict_fields = scalar(keys %$dict_href);
    
    my @column_headers = $self->get_column_headers();
        
    
    my $delim = $self->get_delim();

    my @out_fields;
    for my $column_header (@column_headers) {
        my $field = $dict_href->{$column_header};
        unless (defined $field) {
            confess "Error, missing value for required column field: $column_header";
        }
        if ($field =~ /$delim/) {
            # don't allow any delimiters to contaminate the field value, otherwise it'll introduce offsets.
            $field =~ s/$delim/ /g;
        }
        # also avoid newlines, which will also break the output formatting.
        if ($field =~ /\n/) {
            $field =~ s/\n/ /g;
        }
        
        push (@out_fields, $field);
    }

    my $outline = join("\t", @out_fields);

    my $ofh = $self->get_fh();

    print $ofh "$outline\n";
    
    return;
}


1; #EOM

