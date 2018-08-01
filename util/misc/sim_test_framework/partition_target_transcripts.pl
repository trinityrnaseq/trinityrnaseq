#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$ENV{TRINITY_HOME}/PerlLib/");
use Fasta_reader;
use Cwd;
use Data::Dumper;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use List::Util qw (shuffle);

my $help_flag;

my $ref_trans_fa;
my $MIN_REFSEQ_LENGTH = 100;
my $OUT_DIR = "Seqs_dir";

my $MAX_ISOFORMS = -1;
my $MIN_ISOFORMS = 2;


my $usage = <<__EOUSAGE__;

################################################################################
#
#  * Required:
#
#  --ref_trans|R <string>                 reference transcriptome
#
#  * Common Opts:
#
#  --by_Gene                              target all isoforms of a gene at once.
#                                            (requires multiple isoforms, ignores single-iso genes)
#
#  --out_dir|O <string>                   output directory name (default: $OUT_DIR)
#
#  * Misc Opts:
#
#  --min_refseq_length <int>            min length for a reference transcript 
#                                       sequence (default: $MIN_REFSEQ_LENGTH)
#
#  if --by_Gene:
#
#    --min_isoforms <int>                 default: $MIN_ISOFORMS
#    --max_isoforms <int>                 max number of isoforms to test (default: $MAX_ISOFORMS)
#
#    --longest_isoform_only               restricts to only single longest isoform per gene.
#
#    --restrict_to_genes <string>         file containing lists of gene accessions to restrict to.
#
############################################################################################



__EOUSAGE__

    ;



my $BY_GENE_FLAG = 0;
my $LONGEST_ISOFORM_ONLY_FLAG = 0;

my $restrict_to_genes_file = "";

&GetOptions ( 'h' => \$help_flag,
                            
              # required
              'ref_trans|R=s' => \$ref_trans_fa,
              
              # optional
              'out_dir|O=s' => \$OUT_DIR,
              'min_refseq_length=i' => \$MIN_REFSEQ_LENGTH,
              
              'by_Gene' => \$BY_GENE_FLAG,
              
              'max_isoforms=i' => \$MAX_ISOFORMS,
              'min_isoforms=i' => \$MIN_ISOFORMS,
              
              'longest_isoform_only' => \$LONGEST_ISOFORM_ONLY_FLAG,
              
              'restrict_to_genes=s' => \$restrict_to_genes_file,
);



if ($help_flag) {
    die $usage;
}

unless ($ref_trans_fa) { 
    die $usage;
}


main: {
    
    my $BASEDIR = cwd();
    
    
    if ($ref_trans_fa =~ /\.gz$/) {
        my $unzipped = $ref_trans_fa;
        $unzipped =~ s/\.gz$//g;
        if (! -s $unzipped) {
            &process_cmd("gunzip -c $ref_trans_fa > $unzipped");
        }

        $ref_trans_fa = $unzipped;
    }
    
    my $fasta_reader = new Fasta_reader($ref_trans_fa);
    my %fasta_seqs = $fasta_reader->retrieve_all_seqs_hash();

    my %reorganized_fasta_seqs = &reorganize_fasta_seqs(\%fasta_seqs, $BY_GENE_FLAG);

    
    my %restricted_genes;
    if ($restrict_to_genes_file) {
        my @gene_ids = `cat $restrict_to_genes_file`;
        chomp @gene_ids;
        %restricted_genes = map { + $_ => 1 } @gene_ids;
    }
    
    my $total_counter = 0;

    my @accs = keys %reorganized_fasta_seqs;

    my %seen;

    foreach my $acc (@accs) {
        
        if (%restricted_genes && ! exists $restricted_genes{$acc}) {
            # skipping, not in the restricted list.
            next;
        }
        $seen{$acc} = 1;
        
        chdir $BASEDIR or die "Error, cannot cd to $BASEDIR";
        
            
        my $seq_entries_aref = $reorganized_fasta_seqs{$acc};
        
        my @min_length_targets;
        foreach my $entry (@$seq_entries_aref) {
            
            my ($trans_acc, $seq) = ($entry->{acc},
                                     $entry->{seq});
            
            if (length($seq) >= $MIN_REFSEQ_LENGTH && $seq !~ /[^GATC]/i) {
                push (@min_length_targets, $entry);
            }
        }
        
        unless (@min_length_targets) { 
            print STDERR "No min length targets to pursue for $acc .... skipping.\n";
            next; 
        }

        unless (-d $OUT_DIR) {
            mkdir($OUT_DIR) or die $!;
        }
                
        @min_length_targets = reverse sort {length($a->{seq}) <=> length($b->{seq}) } @min_length_targets;
        
        my $num_total_targets = scalar(@min_length_targets);
        
        if ($BY_GENE_FLAG) {
            
            if ($LONGEST_ISOFORM_ONLY_FLAG) {
                @min_length_targets = shift @min_length_targets;
            }
            else {
                
                if ($num_total_targets < $MIN_ISOFORMS) {
                    next;
                }
                
                if ($num_total_targets > $MAX_ISOFORMS) {
                    
                    @min_length_targets = @min_length_targets[0..($num_total_targets-1)];
                }
            }
            
        }

        &prep_seqs($acc, \@min_length_targets);
        
        $total_counter++;
        if ($total_counter % 100 == 0) {
            print STDERR "\n[$total_counter]\n";
        }
    }
    

    if (%restricted_genes) {
        # ensure we got them all
        for my $seen_acc (keys %seen) {
            if (exists $restricted_genes{$seen_acc}) {
                delete $restricted_genes{$seen_acc};
            }
        }

        if (%restricted_genes) {
            die "Error, missing entries for restricted gene list entries: " . Dumper(\%restricted_genes);
        }
        else {
            print STDERR "-all restricted gene entries identified and reported.\n";
        }
    }
    
    print STDERR "\nDone.\n\n";
    
    exit(0);
}



####
sub prep_seqs {
    my ($acc, $entries_aref) = @_;
    
    print STDERR "\r-processing $acc   ";

    my $num_entries = scalar(@$entries_aref);

    my $basedir = cwd();
    
    my $dir_tok = $acc;
    $dir_tok =~ s/\W/_/g;
    
    my $workdir = "$OUT_DIR/$dir_tok";
    
    unless (-d $workdir) {
        mkdir $workdir or die "Error, cannot mkdir $workdir";
    }

    my $refseqs_fa = "$workdir/refseqs.fa";
    
    # write ref fasta seqs.
    if (! -s "$refseqs_fa") {
        open (my $ofh, ">$refseqs_fa") or die $!;
        foreach my $entry (@$entries_aref) {
            my ($entry_acc, $seq) = ($entry->{acc},
                                     $entry->{seq});
            
            print $ofh ">$entry_acc\n$seq\n";
        }
        close $ofh;
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
sub reorganize_fasta_seqs {
    my ($fasta_seqs_href, $by_gene_flag) = @_;

    my %reorg_fasta;

    foreach my $acc (sort keys %$fasta_seqs_href) {
        
        my $seq = uc $fasta_seqs_href->{$acc};

        my $key = $acc;
        if ($by_gene_flag) {
            if ($acc =~ /^([^;]+);([^;]+)$/) {
                my $trans = $1;
                my $gene = $2;
                $key = $gene;
            }
            elsif ($acc =~ /^([^\|]+)\|([^\|]+)$/) {
                my $gene = $1;
                my $trans = $2;
                $key = $gene;
            }
            else {
                confess "Error, no gene ID extracted from $acc ";
            }
            
        }
        
        push (@{$reorg_fasta{$key}}, { acc => $acc, 
                                       seq => $seq}
            );
    }
    
    return(%reorg_fasta);
    
}

