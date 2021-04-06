#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Fasta_reader;
use Nuc_translator;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Cwd;

my $usage = <<__EOUSAGE__;


##################################################################

$0

##################################################################
#
# Required:
#
#  --transcripts <string>       file containing target transcripts in fasta format
#
# Optional:
#
#  --read_length <int>          default: 76
#
#  --frag_length <int>             default: 300
#
#  --out_prefix <string>        default: 'reads'
#
#  --depth_of_cov <int>         default: 100  (100x targeted base coverage)
#
####
#
#  following wgsim options are pass-through:
#
# Options: 
#    -e FLOAT      base error rate [0.020]
#    -s INT        standard deviation [50]
#    -r FLOAT      rate of mutations [0.0010]
#    -R FLOAT      fraction of indels [0.15]
#    -X FLOAT      probability an indel is extended [0.30]
#    -S INT        seed for random generator [-1]
#    -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [0.05]
#    -h            haplotype mode
#    -Z INT        strand specific mode: 1=FR, 2=RF
#    -D            debug mode... highly verbose
#
#
############################################################################################



__EOUSAGE__

    ;



my $require_proper_pairs_flag = 0;

my $transcripts;
my $read_length = 76;
my $frag_length = 300;
my $help_flag;
my $out_prefix = "reads";
my $depth_of_cov = 100;

&GetOptions ( 'help' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'read_length=i' => \$read_length,
              'frag_length=i' => \$frag_length,
              'out_prefix=s' => \$out_prefix,
              'depth_of_cov=i' => \$depth_of_cov,
              
              
    );



if ($help_flag) {
    die $usage;
}

unless ($transcripts) { 
    die $usage;
}

unless ($out_prefix =~ /^\//) {
    $out_prefix = cwd() . "/$out_prefix";
}

main: {
    
    my $number_reads = &estimate_total_read_count($transcripts, $depth_of_cov, $read_length);
    
    my $cmd = "wgsim-trans -N $number_reads -1 $read_length -2 $read_length "
        . " -d $frag_length "
        . " @ARGV > $out_prefix.log"; # pass-through to wgsim
    
    
    my $token = "wgsim_R${read_length}_F${frag_length}_D${depth_of_cov}";
    if (grep { "-Z" } @ARGV) {
        $token .= "_SS";
    }
    
    
    
    my $left_prefix = "$out_prefix.$token.left";
    my $right_prefix = "$out_prefix.$token.right";
    
    $cmd .= " $transcripts $left_prefix.fq $right_prefix.fq";
    
    &process_cmd($cmd);

    # convert to fasta format
    &process_cmd("$FindBin::Bin/../support_scripts/fastQ_to_fastA.pl -I $left_prefix.fq > $left_prefix.fa");

    &process_cmd("$FindBin::Bin/../support_scripts/fastQ_to_fastA.pl -I $right_prefix.fq > $right_prefix.fa");

    #unlink("$left_prefix.fq", "$right_prefix.fq");

    open(my $ofh, ">$out_prefix.$token.info") or die "Error, cannot write to file: $out_prefix.$token.info";
    print $ofh join("\t", $transcripts, "$left_prefix.fa", "$right_prefix.fa");
    close $ofh;
    
    
    exit(0);

}


####
sub estimate_total_read_count {
    my ($transcripts_fasta_file, $depth_of_cov, $read_length) = @_;

    my $sum_seq_length = 0;
    my $fasta_reader = new Fasta_reader($transcripts);
    while (my $seq_obj = $fasta_reader->next()) {
        my $sequence = $seq_obj->get_sequence();
        $sum_seq_length += length($sequence);
    }
    
    # DOC = num_reads * read_length * 2 / seq_length

    # so, 

    # num_reads = DOC * seq_length / 2 / read_length

    my $num_reads = $depth_of_cov * $sum_seq_length / 2 / $read_length;

    $num_reads = int($num_reads + 0.5);
    
    return($num_reads);
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

