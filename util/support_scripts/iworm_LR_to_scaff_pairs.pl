#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

use Carp;
use Getopt::Long qw(:config posix_default no_ignore_case bundling pass_through);
use Data::Dumper;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Fasta_reader;

my $MIN_CONTAINMENT = 96;
my $KMER_LENGTH = 25;

my $usage = <<__EOUSAGE__;

##############################################################################
#
# Required:
#
#   --blast_outfmt6_wlen <string>   LR_blastn.outfmt6.wLen file
#
#   --iworm_fasta <string>          inchworm.K25.L25.fa
#
#   Optional:
#
#   --min_containment <int>         default: $MIN_CONTAINMENT
#
#   --kmer_len <int>                default: $KMER_LENGTH
#
#   --debug
#
###############################################################################


__EOUSAGE__

    ;

my $help_flag;

my $DEBUG;
my $blast_outfile;
my $iworm_fasta;

&GetOptions ( 'h' => \$help_flag,
              'blast_outfmt6_wlen=s' => \$blast_outfile,
              'iworm_fasta=s' => \$iworm_fasta,
              'min_containment=i' => \$MIN_CONTAINMENT,
              'kmer_len=i' => \$KMER_LENGTH,
              'debug' => \$DEBUG,
    );

if ($help_flag) {
    die $usage;
}

unless ($blast_outfile && $iworm_fasta) {
    die $usage;
}


my %ACC_TO_KMERS_CACHE;
my %TRANS_SEQS;

main: {

    my %LR_to_iworm_hits;

    {
        open (my $fh, $blast_outfile) or die "Error, cannot open file $blast_outfile";
        while (<$fh>) {
            #print;
            chomp;
            my @x = split(/\t/);
            my ($iworm_acc, $LR_acc, 
                $iworm_end5, $iworm_end3,
                $LR_end5, $LR_end3,
                $bitscore,
                $iworm_containment,
                ) = ($x[0], $x[1], 
                     $x[6], $x[7], 
                     $x[8], $x[9],
                     $x[11],
                     $x[13],
                    );
            
            my ($LR_lend, $LR_rend) = sort {$a<=>$b} ($LR_end5, $LR_end3);

            push (@{$LR_to_iworm_hits{$LR_acc}}, { iworm_acc => $iworm_acc,
                                                   lend => $LR_lend,
                                                   rend => $LR_rend,
                                                   iworm_containment => $iworm_containment,
                  });
            
            
            
        }
        close $fh;
    }

    
    my $fasta_reader = new Fasta_reader($iworm_fasta);
    %TRANS_SEQS = $fasta_reader->retrieve_all_seqs_hash();
    

    my %iworm_pairs;
    
    foreach my $LR (keys %LR_to_iworm_hits) {

        my @iworm_hits = sort {$a->{lend}<=>$b->{lend}
                               ||
                                   $a->{rend} <=> $b->{rend}} @{$LR_to_iworm_hits{$LR}};


        if ($DEBUG) {
            # simple report of iworm hits:
            print STDERR "\n// Iworm matches for $LR\n";
            foreach my $iworm_hit (@iworm_hits) {
                print STDERR join("\t", $iworm_hit->{iworm_acc}, $iworm_hit->{lend}, 
                                  $iworm_hit->{rend}, $iworm_hit->{iworm_containment}) . "\n";
            }
        }


        my %full_containments;
        
        for (my $i = 0; $i < $#iworm_hits; $i++) {

            my $iworm_i = $iworm_hits[$i];

            for (my $j = $i + 1; $j <= $#iworm_hits; $j++) {

                my $iworm_j = $iworm_hits[$j];

                if ($iworm_i->{iworm_acc} eq $iworm_j->{iworm_acc}) { next; }

                if ($iworm_j->{lend} > $iworm_i->{rend}) {
                    last; 
                }

                unless ($iworm_i->{iworm_containment} >= $MIN_CONTAINMENT || $iworm_j->{iworm_containment} >= $MIN_CONTAINMENT) {
                    next; 
                }
            
                if ($iworm_i->{iworm_containment} >= $MIN_CONTAINMENT) {
                    $full_containments{ $iworm_i->{iworm_acc} } = 1;
                }
                if ($iworm_j->{iworm_containment} >= $MIN_CONTAINMENT) {
                    $full_containments{ $iworm_j->{iworm_acc} } = 1;
                }
                
                my $overlap_len = &get_overlap_len($iworm_i, $iworm_j);
                
                my $share_kmer_overlap_flag = &share_kmer($iworm_i->{iworm_acc}, $iworm_j->{iworm_acc}); 
                
                if ($DEBUG) {
                    print STDERR join("\t", $iworm_i->{iworm_acc}, $iworm_j->{iworm_acc}, 
                                      "overlap_len:$overlap_len", "share_kmer_overlap:$share_kmer_overlap_flag")  . "\n";
                }
                
                if ($overlap_len && $share_kmer_overlap_flag) {
                    
                    my $pair_token = join("$;", sort ($iworm_i->{iworm_acc}, $iworm_j->{iworm_acc}));
                    $iworm_pairs{$pair_token} = &get_min_iworm_cov($iworm_i->{iworm_acc}, $iworm_j->{iworm_acc});
                    
                }
            }
        } # end of i-vs-j comparisons.

        my @fully_contained_iworms = keys %full_containments;
        if (scalar(@fully_contained_iworms) > 1) {
            # make all pairwise links
            for (my $i = 0; $i < $#fully_contained_iworms; $i++) {
                my $iworm_i_acc = $fully_contained_iworms[$i];
                for (my $j = $i + 1; $j <= $#fully_contained_iworms; $j++) {
                    my $iworm_j_acc = $fully_contained_iworms[$j];
                    my $pair_token = join("$;", sort ($iworm_i_acc, $iworm_j_acc));
                    $iworm_pairs{$pair_token} = &get_min_iworm_cov($iworm_i_acc, $iworm_j_acc);
                }
            }
        }
    }
    

    ## get fasta index values for the iworm contigs:
    my %iworm_to_fasta_index;
    {
        my $index = 0;
        open (my $fh, $iworm_fasta) or die "Error, cannot open file $iworm_fasta";
        while (<$fh>) {
            chomp;
            if (/^>(\S+)/) {
                my $acc = $1;
                $iworm_to_fasta_index{$acc} = $index;
                $index++;
            }
        }
        close $fh;
    }

    
    foreach my $pair (keys %iworm_pairs) {

        my ($iworm_acc_A, $iworm_acc_B) = split(/$;/, $pair);
        
        # a179572;19
        
        # assign glue as the lower of the coverage values:

        my ($coreA, $covA) = split(/;/, $iworm_acc_A);
        my ($coreB, $covB) = split(/;/, $iworm_acc_B);
        
        my $glue = min($covA, $covB);
        
        # output record
        my $iworm_index_A = $iworm_to_fasta_index{$iworm_acc_A};
        
        my $iworm_index_B = $iworm_to_fasta_index{$iworm_acc_B};

                
        print join("\t", 
                   $iworm_acc_A, $iworm_index_A, 
                   $iworm_acc_B, $iworm_index_B,
                   $glue) . "\n";
    }



    exit(0);
    
}


####
sub get_min_iworm_cov {
    my ($acc_i, $acc_j) = @_;

    my ($pref_i, $cov_i) = split(/;/, $acc_i);
    
    my ($pref_j, $cov_j) = split(/;/, $acc_j);

    if ($cov_i < $cov_j) {
        return($cov_i);
    }
    else {
        return($cov_j);
    }
}



####
sub get_overlap_len {
    my ($iworm_i, $iworm_j) = @_;

    ## require:
    ##
    ##   -----------------   iworm_i
    ##             -------------------   iworm_j
    ##
    ## ----------------------------------  LR
    ##
    ##  overlap by at least K-1
    ##  staggered overlaps along LR
    ##  iworm_i and iworm_j alignment lengths >= 2 * (K-1)

    
    my $iworm_i_align_len = $iworm_i->{rend} - $iworm_i->{lend} + 1;
    my $iworm_j_align_len = $iworm_j->{rend} - $iworm_j->{lend} + 1;

    
    if ($iworm_i->{rend} > $iworm_j->{rend}) {
        return(0);
    }
    
    
    my $overlap_len = $iworm_i->{rend} - $iworm_j->{lend} + 1;
    
    return($overlap_len);
}


####
sub share_kmer {
    my ($accA, $accB) = @_;

    my $kmers_href_A = &get_kmers($accA);
        
    my $kmers_href_B = &get_kmers($accB);

    foreach my $kmer (keys %$kmers_href_A) {
        if (exists $kmers_href_B->{$kmer}) {
            return(1); # YES
        }
    }


    return(0); # NO
    
}

####
sub get_kmers {
    my ($acc) = @_;
    
    if (my $href = $ACC_TO_KMERS_CACHE{$acc}) {
        return($href);
    }
    else {
        my $sequence = uc $TRANS_SEQS{$acc};
        
        my %kmers;
        for (my $i = 0; $i <= length($sequence) - ($KMER_LENGTH-1); $i++) {

            my $kmer = substr($sequence, $i, $KMER_LENGTH-1);
            $kmers{$kmer} = 1;
        }
        $ACC_TO_KMERS_CACHE{$acc} = \%kmers;
        
        return($ACC_TO_KMERS_CACHE{$acc});
    }

}
