#!/usr/bin/env perl

use strict;

my $SEE = 0;

###############
# blat format: # Q=cDNA T=genomic
################

#  0: match
#  1: mis-match
#  2: rep. match
#  3: N's
#  4: Q gap count
#  5: Q gap bases
#  6: T gap count
#  7: T gap bases
#  8: strand
#  9: Q name
# 10: Q size
# 11: Q start
# 12: Q end
# 13: T name
# 14: T size
# 15: T start
# 16: T end
# 17: block count
# 18: block Sizes
# 19: Q starts
# 20: T starts
# 21: Q seqs (pslx format)
# 22: T seqs (pslx format)

#############################################################
## This script filters blat output and extracts only       ##
## the top scoring alignment chain for each accession.     ##
#############################################################

my $line_num = 0;
my %data;

my $filename = $ARGV[0] or die "\n\nusage: $0 outputfile.psl [num_top_hits=1]\n\n";
my $num_top_hits = $ARGV[1] || 1;

## First Pass, assign scores to entries associated with accessions.
my $scores_file = "$filename.scores";
{
    
    open (my $ofh, ">$scores_file") or "Error, cannot write to file $scores_file";
    open (FILE, "$filename");
    my $line_num = 0;
    while (<FILE>) {
        $line_num++;
        my @x = split (/\t/);
        unless ($x[0] =~ /^\d/) {next;}
        my ($matches, $mismatches, $q_gap_num, $t_gap_num, $q_insert, $t_insert) = ($x[0], $x[1], $x[4], $x[6], $x[5], abs($x[7]));
        
        my $score = &calculate_score($matches, $mismatches, $q_gap_num, $t_gap_num, $q_insert, $t_insert);
        
        $score = sprintf("%.2f", $score);
       
        my $accession = $x[9];
        
        print $ofh join("\t", $accession, $score, $line_num) . "\n";


        #if ($line_num > 10000) { last; }
        
    }
    
    close FILE;
    close $ofh;
}

my $sort_by_score_file = "$scores_file.sort_by_score";
{
    ## sort the file by accession, score desc
    my $cmd = "sort -k1,1 -k2,2nr $scores_file > $sort_by_score_file";
    &process_cmd($cmd);
}


my $top_line_file = "$sort_by_score_file.top";
{
    open (my $ofh, ">$top_line_file") or die "Error, cannot write to $top_line_file";
    
    ## Identify each highest scoring match
    my $prev_acc = "";
    my $count_reported = 0;
    open (my $fh, $sort_by_score_file) or die $!;
    while (<$fh>) {
        my $line = $_;
        chomp;
        my ($acc, $score, $line_no) = split(/\t/);
        if ($acc ne $prev_acc) {
            ## highest score comes first.
            

            $prev_acc = $acc;
            $count_reported = 0;
        }
        
        if ($count_reported < $num_top_hits) {
            print $ofh $line;
        }

        $count_reported++;
        
        
    }
    close $fh;
    close $ofh;
}
          


my $sort_by_line_no_file = "$top_line_file.sort_by_line_no";
{
    ## now sort by line number
    my $cmd = "sort -k3,3n $top_line_file > $sort_by_line_no_file";
    &process_cmd($cmd);
    
}


open (my $lines_fh, $sort_by_line_no_file) or die $!;
my $line = <$lines_fh>;
chomp $line;
my ($acc, $score, $line_no) = split(/\t/, $line);


open (FILE, "$filename");
$line_num = 0;
while (<FILE>) {
    my $outline = $_;
    $line_num++;
    if ($line_num == $line_no) {
        print $outline;
        
        ## set next search line no
        my $line = <$lines_fh>;
        if (! $line) {
            # no more to search 
            last;
        }
        ($acc, $score, $line_no) = split(/\t/, $line);
        
    }
    elsif ($line_num >$line_no) {
        die "Error, $line_num exceeded $line_no";
    }
}
close FILE;

## cleanup
unlink($scores_file, $sort_by_score_file, $top_line_file, $sort_by_line_no_file);

exit(0);


####
sub calculate_score { 
    my ($match, $mismatch, $q_gap_num, $t_gap_num, $q_insert, $t_insert) = @_;
    
    ## JKent's code in pslFilter:
    # score = (psl->match + psl->repMatch)*reward - psl->misMatch*cost 
    #    - (psl->qNumInsert + psl->tNumInsert + 1) * gapOpenCost
    #	- log(psl->qBaseInsert + psl->tBaseInsert + 1) * gapSizeLogMod;
    
    #use jkent's default score parameters:
    my $reward = 1;
    my $cost = 1;
    my $gapOpenCost = 4;
    my $gapSizeLogMod = 1;
    
    return ( ($match * $reward) - ($mismatch * $cost) - 
             ( ($q_gap_num + $t_gap_num) * $gapOpenCost) - 
             ( log ($q_insert + $t_insert + 1) * $gapSizeLogMod) );
    
} 



####
sub process_cmd {
    my ($cmd) = @_;

    print STDERR "CMD: $cmd\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Error, cmd: $cmd died with ret $ret";
    }
    
    return;
}
