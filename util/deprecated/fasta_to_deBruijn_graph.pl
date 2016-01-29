#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../PerlLib", "$FindBin::RealBin/../PerlLib/KmerGraphLib");
use Nuc_translator;
use KmerGraph;
use Fasta_reader;

use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

#######################################################################
#
# Required:
#
#  --fasta <string>     fasta file (or inchworm bundle file)
#  -K <int>             kmer size  (overlaps by K-1)
#  -C <int>             component identifier
#
# Optional:
#
#  --SS                 indicates strand-specific
#
#######################################################################

__EOUSAGE__

    ;


my ($help_flag, $fasta_file, $kmer_size, $strand_specific_mode);
my $component_id;


&GetOptions ('h' => \$help_flag,
             'fasta=s' => \$fasta_file,
             'K=i' => \$kmer_size,
             'C=i' => \$component_id,
             'SS' => \$strand_specific_mode,
    );

if ($help_flag) {
    die $usage;
}

unless ($fasta_file && $kmer_size && defined($component_id)) { 
    die $usage;
}

main: {

    my $kmer_graph = new KmerGraph($kmer_size);

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {
        
        my $acc = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        
        my @sequences = split(/X/, $sequence); # in case of iworm bundles
        
        foreach my $seq_string (@sequences) {

            $kmer_graph->add_sequence_to_graph($acc, $seq_string, 1, 'black');

            unless ($strand_specific_mode) {
                my $rc_seq = &reverse_complement($seq_string);
                $kmer_graph->add_sequence_to_graph($acc, $rc_seq, 1, 'black');
            }
        }

    }

    # print $kmer_graph->toGraphViz();

    print "Component $component_id\n";
    
    my @nodes = $kmer_graph->get_all_nodes();
    
    


    my @root_nodes = $kmer_graph->get_root_nodes();
    unless (@root_nodes) {
        # must be circular, just choose one node as a root
        @root_nodes = ($nodes[0]);
    }
    
    my @forward_queue;

    my @reports;
    
    foreach my $root_node (@root_nodes) {

        
        my @visited_nodes;
        
        push (@forward_queue, $root_node);
        
        while (@forward_queue) {
            my $node = shift @forward_queue;
            my $node_id = $node->get_ID();
            
            if ($node->{visited}) {
                next;
            }
            
            my $kmer = $node->get_value();
            
            my @prev_node_ids;
            my @prev_nodes;
            if (@prev_nodes = $node->get_all_prev_nodes()) {
                foreach my $prev_node (@prev_nodes) {
                    push (@prev_node_ids, $prev_node->get_ID());
                }
            }
            else {
                @prev_node_ids = (-1);
            }
            
            foreach my $prev_id (@prev_node_ids) {
                my $out = join("\t", $node_id, $prev_id, 1, $kmer, 0);
                push (@reports, [$node_id, $out]);
            }
                
            if (my @next_nodes = $node->get_all_next_nodes()) {
                foreach my $next_node (@next_nodes) {
                    unless ($next_node->{visited}) {
                        push (@forward_queue, $next_node);
                    }
                }
            }
            if (@prev_nodes) {
                foreach my $prev_node (@prev_nodes) {
                    unless ($prev_node->{visited}) {
                        push (@forward_queue, $prev_node);
                    }
                }
            }
            $node->{visited} = 1;
            push (@visited_nodes, $node) unless ($strand_specific_mode);
            
            
        }
        
        unless ($strand_specific_mode) {
            ## deactivate the reverse complements
            foreach my $visited_node (@visited_nodes) {
                my $kmer = $visited_node->get_value();
                my $rc_node = $kmer_graph->get_node(&reverse_complement($kmer));
                unless ($rc_node) {
                    die "Error, no node for rc($kmer) ";
                }
                $rc_node->{visited} = 1;
                # print STDERR "-deactivating: " . $rc_node->get_ID() . "\n";
            }
        }
    }
    

    ## Print graph 

    @reports = sort {$a->[0]<=>$b->[0]} @reports;
    foreach my $set (@reports) {
        my ($node_id, $text) = @$set;
        print "$text\n";
    }
    

    exit(0);
}



        
        
        











