#!/usr/bin/env perl

# works for trinity transcripts or trinity genes, Sept 28 2016 bhaas


# lightweight fasta reader capabilities:
package Fasta_reader;

use strict;

sub new {
    my ($packagename, $fastaFile) = @_;

	## note: fastaFile can be a filename or an IO::Handle
	

    my $self = { fastaFile => undef,,
				 fileHandle => undef };

    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    
	if (ref $fastaFile eq 'IO::Handle') {
		$filehandle = $fastaFile;
	}
	else {
		
		open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
		$self->{fastaFile} = $fastaFile;
	}
	
	$self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;
    
	if (defined($next_text_input) && $next_text_input !~ /\w/) {
		## must have been some whitespace at start of fasta file, before first entry.
		## try again:
		$next_text_input = <$filehandle>;
	}
	
	my $seqobj = undef;
    
	if ($next_text_input) {
		$next_text_input =~ s/^>|>$//g; #remove trailing > char.
		$next_text_input =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove cntrl chars
		my ($header, @seqlines) = split (/\n/, $next_text_input);
		my $sequence = join ("", @seqlines);
		$sequence =~ s/\s//g;
		
		$seqobj = Sequence->new($header, $sequence);
    }
    
    $/ = $orig_record_sep; #reset the record separator to original setting.
    
    return ($seqobj); #returns null if not instantiated.
}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}

####
sub retrieve_all_seqs_hash {
	my $self = shift;

	my %acc_to_seq;
	
	while (my $seq_obj = $self->next()) {
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();

		$acc_to_seq{$acc} = $sequence;
	}

	return(%acc_to_seq);
}



##############################################
package Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;
    
    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
        
    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
		 filename => undef };
    bless ($self, $packagename);
    return ($self);
}

####
sub get_accession {
    my $self = shift;
    return ($self->{accession});
}

####
sub get_header {
    my $self = shift;
    return ($self->{header});
}

####
sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}

#### 
sub get_FASTA_format {
    my $self = shift;
    my $header = $self->get_header();
    my $sequence = $self->get_sequence();
    $sequence =~ s/(\S{60})/$1\n/g;
    my $fasta_entry = ">$header\n$sequence\n";
    return ($fasta_entry);
}


####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;

    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});
    
	my $fasta_entry = $self->get_FASTA_format();
	
    my $tempfile;
    if ($filename) {
		$tempfile = $filename;
    } else {
		my $acc = $accession;
		$acc =~ s/\W/_/g;
		$tempfile = "$acc.fasta";
    }
    
    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP $fasta_entry;
    close TMP;
    return ($tempfile);
}

package main;

my $usage = "usage: $0 acc.list.txt file.fasta\n\n";

my $acc_list = $ARGV[0] or die $usage;
my $pep = $ARGV[1] or die $usage;

main: {
    my $fasta_reader = new Fasta_reader($pep);
    
    my $acc_text = `cat $acc_list`;
    
    my %accs;

    while ($acc_text =~ /(\S+)/g) {
        my $acc = $1;
        
        $accs{$acc} = 1;
    }

    my %seen;
    while (my $seq_obj = $fasta_reader->next()) {

        my $acc = $seq_obj->get_accession();

        my $gene_id = $acc;
        $gene_id =~ s/_i\d+$//;
        
        if ($accs{$acc} || $accs{$gene_id}) {
            print $seq_obj->get_FASTA_format();
            $seen{$acc} = 1 if $accs{$acc};
            $seen{$gene_id} = 1 if $accs{$gene_id};
            
        }
    }

    # remove seen entries
    foreach my $seen_acc (keys %seen) {
        delete $accs{$seen_acc} if exists $accs{$seen_acc};
    }

    
    if (%accs) {
        print STDERR "Error, could not locate entries for: " . join(", ", keys %accs) . "\n";
        exit(1);
    }
    
    exit(0);
}


        
        




