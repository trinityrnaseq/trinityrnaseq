#!/usr/bin/env perl

=pod

=head1 NAME

  Monitor a PBS job progress

=head1 USAGE

  <job identifier> [options]

 Options:
  one of 
   -top         use top on execution host
   -dump        dump output or error (default)
   -follow      follow output/error
   -tail        tail of output/error (the last 10 lines)
   -head        head of output/error (the first 10 lines)

 and also 
   -e|error   Show stderr instead of stdout
   -s|spool   Location of spool directory (defaults to /var/spool/PBS/spool)

=head1 LICENSE

 Released under the MIT License - Alexie Papanicolaou 2012, CSIRO Ecosystem Sciences, alexie@butterflybase.org

=cut


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
my ($head,$follow,$tail,$dump,$show_error,$spool_dir,$top);
GetOptions(
	'top' => \$top,
	'head' => \$head,
	'follow' => \$follow,
	'tail' => \$tail,
	'cat|dump' => \$dump,
	's|spool_dir:s' => \$spool_dir,
	'error' => \$show_error,
);

my $method;
if ($top){
  $method = 'top';
}elsif ($head){
  $method = 'head';
}elsif ($follow){
  $method = 'follow'
}elsif ($tail){
  $method = 'tail';
}else{
  $method = 'cat';
}

my $jobid = shift;

pod2usage unless $jobid;

$spool_dir = $spool_dir ? $spool_dir : '/var/spool/PBS/spool'; # exists on exec host but necessarily on submit host
my $user = $ENV{'USER'};
my $exec = 'ssh ';
pod2usage "No job ID provided\n" unless $jobid;
$jobid=~/^(\w+\[?\d*\]?)/;
$jobid=$1 || pod2usage "Not a valid job ID $jobid\n";

my $pbs_server = `qstat -Bf|grep ^Server`;
chomp($pbs_server);
$pbs_server=~s/^Server: //;
die "No PBS server found. Is it online?\n" unless $pbs_server;

my $node=`qstat -f $jobid|grep exec_host`; 
$node =~/exec_host\s+=\s+(\w+)/;
$node = $1 || "No valid host found. Is $jobid a live job?\n";
my $cmd;

if ($method=~/^c/){
	$cmd = 'cat';
}elsif ($method=~/^ta/){
	$cmd = 'tail';
}elsif ($method=~/^f/){
	$cmd = 'tail -f';
}elsif ($method=~/^h/){
	$cmd = 'head';
}else {
  $cmd = $method;
}
if ($method eq 'top'){
 system("ssh $node -t top");
}else{
  $exec.= " $node $cmd $spool_dir/$jobid.$pbs_server.OU" if !$show_error;
  $exec.= " $node $cmd $spool_dir/$jobid.$pbs_server.ER" if $show_error;
  system($exec);
}
print "\n#\tQPEEK complete for job $jobid. Host was: $node\n";

