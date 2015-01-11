package VCF_parser;

use strict;
use warnings;

use Carp;

sub new {
    my ($packagename) = shift;
    my ($filename) = @_;

    unless ($filename) {
        confess "Error, need filename as parameter";
    }

    my $self = { filename => $filename,
                 fh => undef,
             };

    bless ($self, $packagename);

    $self->_init();
    
    return($self);
}

####
sub _init {
    my ($self) = @_;

    my $filename = $self->{filename};
    open (my $fh, $filename) or confess "Error, cannot open file $filename";
    
    $self->{fh} = $fh;
    
    return;
}


####
sub get_next {
    my $self = shift;
    
    my $fh = $self->{fh};
    my $line = <$fh>;
    while ($line && $line =~ /^\#/) {
        # skip the header lines
        $line = <$fh>;
    }

    if ($line) {
        return(VCF_record->new($line));
    }
    else {
        return(undef);
    }
}




####################################
####################################

package VCF_record;

use strict;
use warnings;

use Carp;

sub new {
    my ($packagename) = shift;
    my ($vcf_line) = @_;

    unless ($vcf_line =~ /\w/) {
        confess "Error, require vcf line of text as parameter";
    }
    
    my $struct = &_parse_vcf_line($vcf_line);
    
    bless($struct, $packagename);

    return($struct);
}

####
sub _parse_vcf_line {
    my ($vcf_line) = @_;
    
    chomp $vcf_line;
        
    my @x = split(/\t/, $vcf_line);
    
    my $acc = $x[0];
    my $pos = $x[1];
    my $ref_base = $x[3];
    my $allele_base = $x[4];
    
    my $tag_info = $x[7];
    
    my %tags;

    foreach my $keyval_pair (split(/;/, $tag_info)) {
        
        if ($keyval_pair =~ /=/) {
            my ($key, $val) = split(/=/, $keyval_pair);
            $tags{$key} = $val;
        }
    }

    my $struct = { line => $vcf_line,
                   
                   acc => $acc,
                   pos => $pos,
                   
                   ref_base => $ref_base,
                   allele_base => $allele_base,
                   tag_info => $tag_info,
                   tags_href => \%tags,
               };

    
    return($struct);
}


####
sub get_accession {
    my ($self) = @_;
    return($self->{acc});
}


####
sub get_position {
    my ($self) = @_;
    return($self->{pos});
}

####
sub get_ref_base {
    my ($self) = @_;
    return($self->{ref_base});
}

####
sub get_allelic_base {
    my ($self) = @_;
    return($self->{allele_base});
}

####
sub get_tag_val {
    my ($self) = shift;
    my ($tagname) = @_;
    
    return($self->{tags_href}->{$tagname});
}

####
sub has_tag {
    my ($self) = shift;
    my ($tagname) = @_;
    
    return(exists $self->{tags_href}->{$tagname});
}




1; #EOM
                   
