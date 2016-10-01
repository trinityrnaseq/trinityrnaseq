package Fastq_reader;

use strict;
use warnings;

sub new {
    my ($packagename, $fastqFile) = @_;

	## note: fastqFile can be a filename or an IO::Handle
	

    my $self = { fastqFile => undef,
				 fileHandle => undef };

    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    
	if (ref $fastqFile eq 'IO::Handle') {
		$filehandle = $fastqFile;
	}
	else {
		if ( $fastqFile =~ /\.gz$/ ) {
		    open ($filehandle, "gunzip -c $fastqFile | ") or die "Error: Couldn't open compressed $fastqFile\n";
        }
        elsif ($fastqFile =~ /\.bz2$/) {
            open ($filehandle, "bunzip2 -c $fastqFile | ") or die "Error, couldn't open compressed $fastqFile $!";
            
        } else {
            open ($filehandle, $fastqFile) or die "Error: Couldn't open $fastqFile\n";
        }
        
		$self->{fastqFile} = $fastqFile;
	}
	
	$self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    
    my $filehandle = $self->{fileHandle};
    my $next_text_input = "";
    
    if (! eof($filehandle)) {
        for (1..4) {
            $next_text_input .= <$filehandle>;
        }
    }
    
	my $read_obj = undef;
    
	if ($next_text_input) {
        
		
		$read_obj = Fastq_record->new($next_text_input);
    }
        
    return ($read_obj); #returns null if not instantiated.

}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}


##############################################
package Fastq_record;

use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my ($text_lines) = @_;
    
    my @split_text = split(/\n/, $text_lines);
    unless (scalar @split_text == 4) {
        confess "Error, fastQ entry doesn't have 4 lines: " . $text_lines;
    }
    
    my ($name_line, $seq_line, $plus, $qual_line) = @split_text;
    
    unless ($name_line =~ /^\@/) { 
        confess "Error, cannot identify first line as read name line: " . $text_lines;
    }
    
    my ($read_name, $rest) = split(/\s+/, $name_line);
    $read_name =~ s/^\@//;
    
    
    my $pair_dir = 0; # assume single

    if ($read_name =~ /^(\S+)\/([12])$/) {
        $read_name = $1;
        $pair_dir = $2;
    }
    elsif (defined($rest) && $rest =~ /^([12]):/) {
        $pair_dir = $1;
    }
    
    
    my $self = { core_read_name => $read_name,
                 pair_dir => $pair_dir, # (0, 1, or 2), with 0 = unpaired.
                 sequence => $seq_line,
                 quals => $qual_line,
                 record => $text_lines,
             };
    

    bless ($self, $packagename);
    return ($self);
}

####
sub get_core_read_name {
    my $self = shift;
    return ($self->{core_read_name});
}

####
sub get_full_read_name {
    my $self = shift;
    
    my $read_name = $self->{core_read_name};
    if ($self->{pair_dir}) {
        return(join("/", $read_name, $self->{pair_dir}));
    }
}

####
sub get_sequence {
    my $self = shift;
    return($self->{sequence});
}

####
sub get_quals {
    my $self = shift;
    return($self->{quals});
}

####
sub get_fastq_record {
    my $self = shift;
    return($self->{record});
}




1; #EOM


