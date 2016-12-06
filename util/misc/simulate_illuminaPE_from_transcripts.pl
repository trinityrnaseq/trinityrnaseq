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
#
# Required:
#
#  --transcripts <string>       file containing target transcripts in fasta format
#
# Optional:
#
#  --read_length <int>          default: 76
#
#  --spacing <int>              default: 1  (simulate read from every (spacing) position)
#
#  --frag_length <int>             default: 300
#  --frag_length_step <int>        default: 100  (only if max_depth > 1)
#
#  --out_prefix <string>        default: 'reads'
#
#  --require_proper_pairs       default(off)
#
#  --include_volcano_spread     default(off)
#
#  --error_rate <float>         default(0), for 1%, set to 0.01
#
#  --max_depth <int>            default(1),  note will be double this if --include_volcano_spread is set.
#
# :: note, generates left.fa and right.fa in FR stranded format
#
#  --make_fastq                 generate fastq instead of fasta (note, uses 'C' for all the qual scores)
#
#################################################################

__EOUSAGE__

    ;



my $require_proper_pairs_flag = 0;

my $transcripts;
my $read_length = 76;
my $spacing = 1;
my $frag_length = 300;
my $frag_length_step = 100;
my $help_flag;
my $out_prefix = "reads";
my $include_volcano_spread = 0;
my $error_rate = 0;
my $MAX_DEPTH = 1;

my $make_fastq_flag = 0;

&GetOptions ( 'h' => \$help_flag,
              'transcripts=s' => \$transcripts,
              'read_length=i' => \$read_length,
              'spacing=i' => \$spacing,
              'frag_length=i' => \$frag_length,
              'frag_length_step=i' => \$frag_length_step,
              'out_prefix=s' => \$out_prefix,
              'require_proper_pairs' => \$require_proper_pairs_flag,
              'include_volcano_spread' => \$include_volcano_spread,
              'error_rate=f' => \$error_rate,
              'max_depth=i' => \$MAX_DEPTH,
              'make_fastq' => \$make_fastq_flag,
    
    );


if ($help_flag) {
    die $usage;
}

unless ($transcripts) { 
    die $usage;
}

if ($error_rate > 0.25) {
    die "Error, error rate is set to $error_rate, which exceeds a max of 0.25 ";
}


main: {
    
    my $fasta_reader = new Fasta_reader($transcripts);
    print STDERR "-parsing incoming $transcripts...";
    my %read_seqs = $fasta_reader->retrieve_all_seqs_hash();
    print STDERR "done.\n";
    
    my $num_trans = scalar (keys %read_seqs);
    my $counter = 0;
    
    my $ext = ($make_fastq_flag) ? "fq" : "fa";

    unless ($out_prefix =~ /^\//) {
        $out_prefix = cwd() . "/$out_prefix";
    }
    
    $out_prefix = "$out_prefix.simPE_R${read_length}_F${frag_length}_FR";
    
    open (my $left_ofh, ">$out_prefix.left.$ext") or die $!;
    open (my $right_ofh, ">$out_prefix.right.$ext") or die $!;
    
    {
        # write info file
        open (my $ofh, ">$out_prefix.info") or die "Error, cannot write to $out_prefix.info";
        print $ofh join("\t", $transcripts, "$out_prefix.left.$ext", "$out_prefix.right.$ext") . "\n";
        close $ofh;
    }
    
    my $FRAG_LENGTH = $frag_length;

    
    foreach my $read_acc (keys %read_seqs) {

        my $seq = uc $read_seqs{$read_acc};

        $frag_length = $FRAG_LENGTH; ## reinit
        
        for my $depth (1..$MAX_DEPTH) {
                        
            $counter++;
            print STDERR "\r[" . sprintf("%.2f%%  = $counter/$num_trans]     ", $counter/$num_trans*100);
                        
            ## uniform dist
            for (my $i = 0; $i <= length($seq); $i+=$spacing) {
                
                my $left_read_seq = "";
                my $right_read_seq = "";
                my $ill_acc = $read_acc . "_Ap$i-D$depth-F$frag_length-$counter";
                
                my $left_start = $i;
                if ($left_start >= 0) {
                    $left_read_seq = substr($seq, $left_start, $read_length);
                }
                
                my $right_start = $i + $frag_length - $read_length + 1;
                if ($right_start + $read_length  <= length($seq)) {
                    $right_read_seq = substr($seq, $right_start, $read_length);
                }
                
                if ($require_proper_pairs_flag && ! ($left_read_seq && $right_read_seq)) { next; }
                
                
                if ($left_read_seq) {
                    if ($error_rate > 0) {
                        $left_read_seq = &introduce_errors($left_read_seq, $error_rate);
                    }
                    
                    &write_seq($left_ofh, "$ill_acc/1", $left_read_seq);
                }
                if ($right_read_seq) {
                    $right_read_seq = &reverse_complement($right_read_seq);
                    
                    if ($error_rate > 0) {
                        $right_read_seq = &introduce_errors($right_read_seq, $error_rate);
                    }
                    
                    &write_seq($right_ofh, "$ill_acc/2", $right_read_seq);
                }


                
            }
            $frag_length += $frag_length_step;
        } # end depth
            
            
      volcano_spread:
        if ($include_volcano_spread) {
            
            ## volcano spread
            for (my $i = 0; $i <= length($seq)/2; $i+=$spacing) {
                
                my $left_read_seq = "";
                my $right_read_seq = "";
                my $ill_acc = $read_acc . "_Bp$i-F$frag_length";
                
                
                my $left_start = $i;
                if ($left_start >= 0) {
                        $left_read_seq = substr($seq, $left_start, $read_length);
                }
                
                my $right_start = length($seq)-$read_length -$i + 1;
                
                if ($left_start + $read_length >= $right_start) { next; } ## don't overlap them.
                
                if ($right_start + $read_length  <= length($seq)) {
                    $right_read_seq = substr($seq, $right_start, $read_length);
                }
                
                
                if ($require_proper_pairs_flag && ! ($left_read_seq && $right_read_seq)) { next; }
                
                if ($left_read_seq) {
                    
                    if ($error_rate > 0) {
                        $left_read_seq = &introduce_errors($left_read_seq, $error_rate);
                    }
                    
                    &write_seq($left_ofh, "$ill_acc/1", $left_read_seq);
                    
                }
                if ($right_read_seq) {
                    
                    $right_read_seq = &reverse_complement($right_read_seq);
                    
                    if ($error_rate > 0) {
                        $right_read_seq = &introduce_errors($right_read_seq, $error_rate);
                    }
                    
                    &write_seq($right_ofh, "$ill_acc/2", $right_read_seq);
                    
                }
                
                
            }
        }
        
    }
    
    close $left_ofh;
    close $right_ofh;
    
    print STDERR "\nDone.\n";
    
    exit(0);
}




####
sub write_seq {
    my ($ofh, $acc, $seq) = @_;

    if ($make_fastq_flag) {
        print $ofh join("\n", "\@$acc", $seq, "+", 'C' x length($seq)) . "\n";
    }
    else {
        # fasta
        print $ofh join("\n", ">$acc", $seq) . "\n";
    }
    return;
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


####
sub capture_kmer_cov_text {
    my ($kmer_cov_file) = @_;
    
    my %kmer_cov;
    
    my $acc = "";
    open (my $fh, $kmer_cov_file) or die "Error, cannt open file $kmer_cov_file";
    while (<$fh>) {
        chomp;
        if (/>(\S+)/) {
            $acc = $1;
        }
        else {
            $kmer_cov{$acc} .= " $_";
        }
    }
    close $fh;

    return(%kmer_cov);
}


####
sub avg {
    my (@vals) = @_;

    if (scalar(@vals) == 1) {
        return($vals[0]);
    }
    

    my $sum = 0;
    foreach my $val (@vals) {
        $sum += $val;
    }
    
    my $avg = $sum / scalar(@vals);


    return(int($avg+0.5));
}

####
sub introduce_errors {
    my ($sequence, $rate) = @_;

    my @seq_chars = split(//, uc $sequence);
    
    my $num_errors = int(length($sequence) * $rate + 0.5);

    if ($num_errors > length($sequence)) {
	confess "ERROR, error rate $rate is yielding more errors than the length of the sequence itself";
    }

    my @mut_chars = qw(G A T C);

    for (1..$num_errors) {

	my $rand_pos = int(rand(length($sequence)));
	
	my $seq_char = $seq_chars[$rand_pos];

	my @possible_mut_chars = grep { $_ ne $seq_char } @mut_chars;

	my $mut_char = $possible_mut_chars[ int(rand(3)) ];
	if ($mut_char eq $seq_char) {
	    die "Error, mut_char == seq_char, shouldn't happen";
	}
	
	$seq_chars[$rand_pos] = $mut_char;
	
    }

    my $mutated_seq = join("", @seq_chars);

    return($mutated_seq);
}


