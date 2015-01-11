#!/usr/bin/perl

use warnings;
use strict;

# see http://perldoc.perl.org/Pod/Usage.html
use Pod::Usage; ## uses pod documentation in usage code

# see http://perldoc.perl.org/Getopt/Long.html
use Getopt::Long qw(:config auto_version auto_help); # for option parsing

# see http://perldoc.perl.org/perlipc.html
use IPC::Open2; # allows safe-ish pipes within perl

=pod

=head1 NAME

readsToComponents.pl -- Uses Bowtie2 to determine the component that
each read maps to

=head1 SYNOPSIS

readsToComponents.pl [options] -i <reads> -f <components>

Additional options :

  -t <int>             : number of CPUs to use
  -o <string>          : output base directory
  -p <int>             : MAPQ minimum value for mapping (default 2)
  -strand              : use single-stranded mapping
  -max_mem_reads <int> : maximum memory to use for sorting reads

=head1 DESCRIPTION

This is a wrapper script (part of the Trinity suite for de-novo
transcriptome assembly) for using Bowtie2 in local mode to determine
the component that each read maps to. It is assumed that a particular
read can only be assigned to a single component, because the greedy
algorithm used for Inchworm enforces that condition on kmers.

The output of this program is a sorted TSV file (component, read,
quality, sequence) containing the reads that map to the
components. Some reads are filtered out in this process:

* MAPQ less than specified lower limit
* reads that are soft-clipped at both ends
* reads containing INDELs relative to the matching component

=cut

my %params = ('outdir' => '.', 'mapq' => '2');

GetOptions(\%params, 'readfile|i=s','compfile|f=s','cpus|t=i',
           'outdir|o=s','strand','mapq|p=i','max_mem_reads=i') or pod2usage();

if(!$params{"readfile"} || !$params{"compfile"}){
  pod2usage("Error: both read file and component file must be specified");
}

# normalise output dir to make sure it ends with '/'
$params{"outdir"} =~ s#/?$#/#;

# /r switch means "return result, but don't change the variable"
my $readExtension = ($params{"readfile"} =~ s/^.*?\.([a-z]+)$/$1/r);
my $bundledBase = ($params{"compfile"} =~ s/\.[a-z]+$//r);
$bundledBase =~ s#^.*/##; # remove directory name from bundled base name

my $readFasta = ($readExtension =~ /fa(sta)?/);

$\ = $/; # automatically add line ending to print statements (i.e. perl -l)

my $bamFileBase = $params{"outdir"}.$bundledBase."_mapped";

warn("-- Generating Bowtie2 index\n");

my @bt2BuildOpts = ('-f', $params{"compfile"},
                    $params{"outdir"}.$bundledBase);

warn("-- bowtie2-build command line: bowtie2-build ".join(" ",@bt2BuildOpts)."\n");

system('bowtie2-build', @bt2BuildOpts) == 0
  or die("Bowtie2 index generation process failed");

warn("-- Mapping reads to components with Bowtie2\n");

# populate program option arrays
# bowtie2: seed size of 25, local search mode

my @bt2Opts = ('-L', '25', '--local',
               '-x', $params{"outdir"}.$bundledBase,
               '-U', $params{"readfile"});
if($readFasta){
  push(@bt2Opts, '-f');
}

# set up bowtie2 run / filter

my $tcOutFileName = $params{"outdir"}. "readsToComponents.out";
open(my $tcOutFile, '>', $tcOutFileName);

warn("-- will write read components to $tcOutFileName\n");


warn("-- Bowtie2 command line: bowtie2 ".join(" ",@bt2Opts)."\n");
# pipe output from bowtie2 run into file handle
my $pidBt2 = open(my $bt2Output, '-|', 'bowtie2', @bt2Opts);

my $totalReads = 0;
my $filteredReads = 0;

while(<$bt2Output>){
  if(!/^@/){ # ignore header lines
    $totalReads++;
    my @F = split(/\t/, $_, 12); # split into SAM fields
    # ignore low mapq, double-ended soft-clipping and INDELs
    if(($F[4] > $params{"mapq"}) && ($F[5] !~ /S.*S/) && ($F[11] =~ /XO:i:0/)){
      $filteredReads++;
      $F[0] = '>'.$F[0];
      $F[2] = substr($F[2],2); # remove s_ from component
      print($tcOutFile join("\t", @F[(2,0,4,9)]));
    }
  }
}

printf(STDERR "-- done (processed %d reads; %d remain after filtering)\n",
      $totalReads, $filteredReads);

close($bt2Output);
close($tcOutFile);
