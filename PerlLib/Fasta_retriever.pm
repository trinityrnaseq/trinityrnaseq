package Fasta_retriever;

use strict;
use warnings;
use Carp;
use threads;
use threads::shared;


my $LOCKVAR :shared;
our $DEBUG = 0;

sub new {
    my ($packagename) = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "Error, need filename as param";
    }

    my $self = { filename => $filename,
                 acc_to_pos_index => undef,
                 fh => undef,
    };
    
    my %acc_to_pos_index :shared;

    $self->{acc_to_pos_index} = \%acc_to_pos_index;
        
    bless ($self, $packagename);

    $self->_init();


    return($self);
}


sub _init {
    my $self = shift;
    
    my $filename = $self->{filename};
    
    # use a samtools faidx index if available
    my $index_file = "$filename.fai";
    if (-s $index_file) {
        open(my $fh, $index_file) or die "Error, cannot open file: $index_file";
        while(<$fh>) {
            chomp;
            my @x = split(/\t/);
            my $acc = $x[0];
            my $file_pos = $x[2];
            $self->{acc_to_pos_index}->{$acc} = $file_pos;
        }
        close $fh;
    }
    else {
        print STDERR "-missing faidx file: $index_file, extracting positions directly.\n";
        print STDERR "-Fasta_retriever:: begin initializing for $filename\n";
    
        open (my $fh, $filename) or die $!;
        $self->{fh} = $fh;
        while (<$fh>) {
            if (/>(\S+)/) {
                my $acc = $1;
                my $file_pos = tell($fh);
                $self->{acc_to_pos_index}->{$acc} = $file_pos;
            }
        }
        print STDERR "-Fasta_retriever:: done initializing for $filename\n";
    }
    
    return;
}

sub refresh_fh {
    my $self = shift;
    
    open (my $fh, $self->{filename}) or die "Error, cannot open file : " . $self->{filename};
    $self->{fh} = $fh;
    
    return $fh;
}


sub get_seq {
    my $self = shift;
    my $acc = shift;

    unless (defined $acc) {
        confess "Error, need acc as param";
    }


    my $seq = "";

    {
        lock $LOCKVAR;
    
        my $file_pos = $self->{acc_to_pos_index}->{$acc} or confess "Error, no seek pos for acc: $acc";
        
        my $fh = $self->refresh_fh();
        seek($fh, $file_pos, 0);

        print STDERR "seeking $acc -> $file_pos\n" if $DEBUG;
        
        while (<$fh>) {
            if (/^>/) {
                print STDERR "   reached $_, stopping\n" if $DEBUG;
                last;
            }
            my $seq_part = $_;
            $seq_part =~ s/\s//g;
            $seq .= $seq_part;
        }
        print STDERR "-done seeking $acc\n\n" if $DEBUG;
        
    }    
    
    return($seq);
}

    

1; #EOM
    
