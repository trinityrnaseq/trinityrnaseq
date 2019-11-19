#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../../PerlLib",  "$FindBin::RealBin/../../PerlLib/KmerGraphLib");
use Fasta_reader;
use ColorGradient;
use Data::Dumper;


my $usage = "usage: $0 graph.seqinfo refseqs.fa\n\n";

my $graph_seqinfo_file = $ARGV[0] or die $usage;
my $refseqs_fa_file = $ARGV[1] or die $usage;

main: {

    my $fasta_reader = new Fasta_reader($refseqs_fa_file);
    my %refseq_fa = $fasta_reader->retrieve_all_seqs_hash();

    my ($node_seq_to_node_id_href, $edges_aref) = &parse_seqinfo($graph_seqinfo_file);
    
    my %refseq_to_nodes_list;
    foreach my $refseq_acc (keys %refseq_fa) {
        my $refseq_seq = $refseq_fa{$refseq_acc};

        foreach my $node_seq (keys %$node_seq_to_node_id_href) {
            my $node_id_info = $node_seq_to_node_id_href->{$node_seq};
            my $idx = &get_node_start_pos($refseq_seq, $node_seq);
            
            if ($idx >= 0) {

                my ($node_id, @rest) = split(/\s+/, $node_id_info);
                
                push (@{$refseq_to_nodes_list{$refseq_acc}}, { start_pos => $idx,
                                                               node_id => $node_id,
                                                               node_info => $node_id_info,
                                                               
                      } );
            }
        }
    }

    my $dot_file = "$graph_seqinfo_file.dot";
    open(my $ofh, ">$dot_file") or die "Error, cannot write to file: $dot_file";
    
    # print Dumper(\%refseq_to_nodes_list);

    ## output graph:
    print $ofh "digraph G {\n" 
        . "    node [width=0.1,height=0.1,fontsize=10];\n"
        . "    edge [fontsize=12];\n"
        . "    margin=1.0;\n"
        . "    rankdir=LR;\n"
        . "    labeljust=l;\n";

        
    
    
    foreach my $node_seq (keys %$node_seq_to_node_id_href) {
        my $node_info_txt = $node_seq_to_node_id_href->{$node_seq};
        print $ofh "    $node_info_txt\n";
    }
    
    ## get a different color for each refseq acc:
    my @colors = &ColorGradient::convert_RGB_hex(&ColorGradient::get_RGB_gradient(scalar keys %refseq_fa));
    
    foreach my $acc (keys %refseq_to_nodes_list) {
        my @structs = @{$refseq_to_nodes_list{$acc}};
        @structs = sort {$a->{start_pos}<=>$b->{start_pos}} @structs;

        my $color = shift @colors;
        
        ## examine each edge
        for (my $i = 0; $i < $#structs; $i++) {
            my $node_id_begin = $structs[$i]->{node_id};
            my $node_id_end = $structs[$i+1]->{node_id};

            print $ofh "    $node_id_begin->$node_id_end [label=\"$acc\", color=\"$color\"];\n";
        }
    }
    
    foreach my $edge (@$edges_aref) {
        print $ofh "    $edge;\n";
    }

    print $ofh "}\n";

    close $ofh;


    exit(system("dot -Tpdf $dot_file > $dot_file.pdf"));
    
    
    
}

####
sub parse_seqinfo {
    my ($seqinfo_file) = @_;

    my %node_seq_to_node_id;
    my @edges;
    
    open(my $fh, $seqinfo_file) or die "Error, cannot open file: $seqinfo_file";
    while(<$fh>) {
        chomp;
        my @x = split(/\t/);
        if (scalar(@x) == 1) {
            my $edge = $x[0];
            if ($edge =~ /(\d+)->(\d+)/) {
                push(@edges, $edge);
            }
            else {
                die "Error, cannot identify $edge as an edge";
            }
        }
        else {
            my ($node_descr, $node_seq) = @x;
            $node_seq_to_node_id{$node_seq} = $node_descr;
        }
    }
    
    close $fh;
    
    return(\%node_seq_to_node_id, \@edges);
}

####
sub get_node_start_pos {
    my ($refseq_seq, $node_seq) = @_;

    my $idx = index($refseq_seq, $node_seq);
    if ($idx >= 0) {
        return($idx);
    }
    else {
        ## try 1st kmer
        my $first_kmer = substr($node_seq, 0, 25);
        $idx = index($refseq_seq, $first_kmer);
        if ($idx) {
            return($idx);
        }
        else {
            # try last kmer
            my $last_kmer = substr($node_seq, -25);
            $idx = index($refseq_seq, $last_kmer);
            return($idx);
        }
    }
}


