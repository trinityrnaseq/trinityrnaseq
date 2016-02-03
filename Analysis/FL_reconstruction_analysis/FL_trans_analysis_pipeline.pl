#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);

use FindBin;
use File::Basename;

$ENV{PATH} .= ":$FindBin::RealBin/../../trinity-plugins/slclust/bin/";

my $help_flag;

my $min_per_id = 99;
my $max_per_gap = 1;
my $min_per_length = 100;


my $usage = <<__EOUSAGE__;

############################################################################################################################
#
# 
#  --target <string>         transcript fasta file (headers must have accessions that look like so:   transcript_id;gene_id
#  --query <string>          Trinity.fa
#  
#  --min_per_id <int>        min percent identity (default: $min_per_id)
#  --max_per_gap <int>        max percent gaps  (default: $max_per_gap)
#  --min_per_length <float>   min percent of reference transcript length to align and still consider FL (default: $min_per_length)
#
#  --forward_orient          strand-specific assemblies, only consider the forward strand.
#
#  --include_tiers           write the top tiers file
#
#  --allow_non_unique_mappings    a single query transcript counts towards all transcripts encapsulated  
#
#  --reuse                   no prompting to reuse any existing blat output files, just use them.
#  --no_reuse                no prompting, and do NOT reuse any existing files.
#
#  --out_prefix              prefix for output files generated.
#
#############################################################################################################################


__EOUSAGE__

    ;


my $target;
my $query;


my $forward_orient = 0;
my $include_tiers = 0;
my $reuse_flag = 0;
my $no_reuse_flag = 0;


my $allow_non_unique_mappings = 0;

my $out_prefix;

&GetOptions ( 'h' => \$help_flag,
              'target=s' => \$target,
              'query=s' => \$query,
              'min_per_id=i' => \$min_per_id,
              'forward_orient' => \$forward_orient,
              'max_per_gap=i' => \$max_per_gap,
              'include_tiers' => \$include_tiers,
              'allow_non_unique_mappings' => \$allow_non_unique_mappings,
              'min_per_length=f' => \$min_per_length,
	      
              'reuse' => \$reuse_flag,
              'no_reuse' => \$no_reuse_flag,
              
              'out_prefix=s' => \$out_prefix,
              );


unless ($target && $query) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

my $util_dir = "$FindBin::RealBin/util";

my $found_all_tools = 1;
my @required_tools = qw (blat slclust);
foreach my $tool (@required_tools) {
    my $path = `which $tool`;
    unless ($path =~ /\w/) {
        print STDERR "Error, cannot locate required tool: $tool\n";
        $found_all_tools = 0;
    }
}
unless ($found_all_tools) {
    die "Error, cannot proceed due to tools missing from PATH setting. ";
}



main: {

    unless ($out_prefix) {
        $out_prefix = basename($query);
    }
    
    my $cmd = "$util_dir/blat_full_length_mappings.pl --target $target --query $query "
        . " --min_per_id $min_per_id --min_per_length $min_per_length --max_per_gap $max_per_gap ";
    if ($forward_orient) {
        $cmd .= " --forward_orient ";
    }
    if ($allow_non_unique_mappings) {
        $cmd .= " --allow_non_unique_mappings ";
    }
    if ($reuse_flag) {
        $cmd .= " --reuse ";
    }
    elsif ($no_reuse_flag) {
        $cmd .= " --no_reuse ";
    }
    

    $cmd .= " --out_prefix $out_prefix ";
    
    $cmd .= " > $out_prefix.pslx.maps";


    &process_cmd($cmd);  ## generates:  $query.pslx,  $query.pslx.FL_entries, $query.pslx.maps
    
    if ($include_tiers) {
        ## generate top tiers
        $cmd = "$util_dir/blat_top_tier_genes.pl $out_prefix.pslx > $out_prefix.pslx.tiers";
        &process_cmd($cmd);
    }
    
    ## generate summary statistics:
    $cmd = "$util_dir/blat_map_filter_with_isoforms.pl $out_prefix.pslx.maps | tee $out_prefix.pslx.maps.summary";
    &process_cmd($cmd);
    
    
    exit(0);

}

####
sub process_cmd {
    my ($cmd) = @_;

    print "CMD: $cmd\n";

    my $ret = system($cmd);

    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }

    return;
}


