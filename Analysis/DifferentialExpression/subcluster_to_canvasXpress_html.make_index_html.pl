#!/usr/bin/env perl

use strict;
use warnings;

use CGI;
use Cwd;

my $cgi = new CGI();
print $cgi->start_html(-title => cwd());

my @entries = `wc -l *.matrix`;
chomp @entries;

print "<table>\n";
print "<tr><th>Num entries</th><th>file</th></tr>\n";
foreach my $entry (@entries) {
    $entry =~ s/^\s+//;
    my ($count, $filename) = split(/\s+/, $entry);
    $count--; # header doesn't count
    print "<tr><td>$count</td><td><a href=\'$filename.html\'>$filename</a></td>\n";
}
print "</table>\n";

print $cgi->end_html();

exit(0);

