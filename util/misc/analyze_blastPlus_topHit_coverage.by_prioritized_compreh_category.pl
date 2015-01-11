#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $help_flag;

my $usage = <<__EOUSAGE__;

#################################################################################################
#
# Required:
#
#  --details <string>              compreh_init_build.details filename
#  --blast <string>                blast.outfmt6.w_pct_hit_length filename
#
# Optional:
#
#  --min_pct_detail_total <float>         minimum percent of total assignments to capture for a given detail category (default: 1.0)
#  --min_pct_hit_len <float>     min length of alignment coverage of hit (default: 80)
#  --max_Evalue <float>          default (1e-20);
#
#################################################################################################


__EOUSAGE__

    ;


my $compreh_details_file;
my $blast_file;
my $min_pct_hit_len = 80;
my $min_pct_detail_total = 1.0;
my $max_Evalue = 1e-20;


&GetOptions ( 'h' => \$help_flag,
              'details=s' => \$compreh_details_file,
              'blast=s' => \$blast_file,
              'min_pct_hit_len=f' => \$min_pct_hit_len,
              'min_pct_detail_total=f' => \$min_pct_detail_total,
              'max_Evalue=f' => \$max_Evalue,
              
              );


if ($help_flag) {
    die $usage;
}

unless ($compreh_details_file && $blast_file) {
    die $usage;
}



my %priorities = ( 'pasa' => 1,
                   'InvalidQualityAlignment' => 2,
                   'PoorAlignment' => 3,
                   'TDN' => 4,
                   );

my %acc_to_detail;
{
    open (my $fh, $compreh_details_file) or die $!;
    while (<$fh>) {
        chomp;
        my ($gene, $trans, $detail) = split(/\t/);
        my ($detail_prefix, $rest) = split(/_/, $detail);
        $acc_to_detail{$trans} = $detail_prefix;
    }
    close $fh;
}


my %FL_mappings; # hit acc => detail
my %hit_to_OS;
my %hit_to_query;

{
    open (my $fh, $blast_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        
        chomp;
        my @x = split(/\t/);
        my $pct_hit_len = $x[13];
        if ($pct_hit_len < $min_pct_hit_len) {
            next;
        }
        
        my $evalue = $x[10];
        if ($evalue > $max_Evalue) { next; }
        
        
        my $query = $x[0];
        my $hit = $x[1];
        my $annot = $x[14];
        
        my $query_detail = $acc_to_detail{$query};
        
        
        if ($annot =~ /OS=(\S+)/) {
            my $os = $1;
            if ($hit_to_OS{$hit}) {
                
                ## see if this has higher priority
                my $prev_detail = $FL_mappings{$hit};
                if ($priorities{$prev_detail} > $priorities{$query_detail}) {
                    $FL_mappings{$hit} = $query_detail; # assign the lower value, meaning higher priority here.
                    $hit_to_query{$hit} = $query;
                }
                
            }
            else {
                $hit_to_OS{$hit} = $os;
                $FL_mappings{$hit} = $query_detail;
                $hit_to_query{$hit} = $query;
            }
        }
        else {
            print STDERR "Error, no OS specified for $annot";
        }
    }
    close $fh;
    
}


## reformat data structure to obtain  detail -> OS -> count
my %data;
my %detail_count_totals;
my %os_total_counts;
my %detail_to_query_list;

{
    foreach my $hit (keys %FL_mappings) {
        
        my $detail = $FL_mappings{$hit};
        my $os = $hit_to_OS{$hit};
        
        my $query = $hit_to_query{$hit};
        
        $data{$detail}->{$os}++;
        $detail_count_totals{$detail}++;

        push (@{$detail_to_query_list{$detail}}, { query => $query, 
                                                   os => $os,
                                               });
    }
    
    open (my $ofh, ">prioritized_query_to_os_summary.len$min_pct_hit_len.E$max_Evalue.txt") or die $!;
    
    ## report
    foreach my $detail (sort {$priorities{$a}<=>$priorities{$b}} keys %priorities) {
        
        my $detail_total = $detail_count_totals{$detail};
        
        print "DETAIL: $detail Total: $detail_total\n";
        
        my $os_counts_href = $data{$detail};
        
        foreach my $os (reverse sort {$os_counts_href->{$a}<=>$os_counts_href->{$b}} keys %$os_counts_href) {
            
            my $count = $os_counts_href->{$os};
            
            my $percent = sprintf("%.2f", $count/$detail_total*100);
            
            if ($percent >= $min_pct_detail_total) {

                print "$os\t$count\t$percent%\n";
                
                $os_total_counts{$os}+= $count;
            }
                


        }
        
        print "\n";
        
        my @entries = @{$detail_to_query_list{$detail}};
        
        foreach my $entry (@entries) {
            
            my $query = $entry->{query};
            my $os = $entry->{os};
            print $ofh join("\t", $detail, $query, $os) . "\n";
        }
        
    }
    close $ofh;
}

## output matrix
my @os_rows = reverse sort {$os_total_counts{$a}<=>$os_total_counts{$b}} keys %os_total_counts;
my @details = sort {$priorities{$a}<=>$priorities{$b}} keys %priorities;


foreach my $type ('counts', 'percentages') {

    print "\n\n";
    print "** Matrix of $type **\n\n";

    print "#\t" . join("\t", @details) . "\n";
    foreach my $os (@os_rows) {
        print "$os";
        foreach my $detail (@details) {
            my $count = $data{$detail}->{$os} || 0;
            
            my $detail_total = $detail_count_totals{$detail};
            my $percent = sprintf("%.2f", $count/$detail_total*100);
            if ($type eq 'counts') {
                print "\t$count";
            }
            else {
                print "\t$percent";
            }
        }
        #my $os_total = $os_total_counts{$os};
        #print "\t\ttotal: $os_total\n";
        print "\n";
    }
    print "\n\n";
}

exit(0);


