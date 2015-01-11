#!/usr/bin/env perl

unless (@params = @ARGV) {
    die "usage: $0 [-t|-s] <number list | range>\n";
}

$delimeter = pop @params;

if ($delimeter eq '-s' || $delimeter eq '-t') {
    if ($delimeter eq '-s') {
	$delimeter = '\s+';
    } elsif ($delimeter eq '-t') {
	$delimeter = '\t';
    }
    pop @ARGV;
} else {
    $delimeter = '\t';
}

if ($ARGV[0] =~ /(^\d+)\-(\d+)/) {
    #print "$1\t$2\n";
    @array = ($1 .. $2);
    
} else {
    @array = @ARGV;
}
    

foreach $entry (@array) {
    $here{$entry} = 1;
}


while (<STDIN>) {
    chomp;
    #my $tab = 0;
    @columns = split (/$delimeter/, $_);
    my $output = "";
    for ($i = 0; $i <= $#columns; $i++) {
		#if ($tab) { print "\t";}
		if ($here{$i}) {
			#print $columns[$i];
			$output .= "$columns[$i]\t";
			#$tab = 1;
		}
    }
	$output =~ s/\s+$//; # remove trailing ws.
    print "$output\n";
}

