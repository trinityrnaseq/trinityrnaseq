#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 chrysalis_dir/ component_no (WELDS|SCAFF|BOTH)\n\n";

my $chrysalis_dir = $ARGV[0] or die $usage;
my $component_no = $ARGV[1];
my $draw_type = $ARGV[2] or die $usage;

unless (defined $component_no) {
    die $usage;
}

unless ($draw_type =~ /^(WELDS|SCAFF|BOTH)$/) {
    die $usage;
}

main: {
    
    my $graph_from_iworm_file = "$chrysalis_dir/GraphFromIwormFasta.out";
    
    my %iworm_accs = &get_iworm_list_from_component($graph_from_iworm_file, $component_no);
    
    ## get the weld links:
    my %welds = &parse_welds($graph_from_iworm_file, \%iworm_accs);
    
    ## get the scaffolding links:
    my $scaffolding_file = "$chrysalis_dir/../iworm_scaffolds.txt";
    my %scaffolds = &parse_scaffolded_pairs($scaffolding_file, \%iworm_accs);
    
    ## write dot file
    print "digraph {\n";
    
    
    my %iworm_to_id;
    ## write node labels
    my $counter = 0;
    foreach my $iworm_acc (keys %iworm_accs) {
        $counter++;
        my $iworm_id = $iworm_to_id{$iworm_acc} = $counter;
        print "\t$iworm_id [label=\"$iworm_acc\"];\n";
    }
 
    if ($draw_type =~ /WELDS|BOTH/) {
        ## add weld links
        foreach my $iworm_acc (keys %welds) {
            my $iworm_id = $iworm_to_id{$iworm_acc};
            foreach my $welded_iworm_acc (keys %{$welds{$iworm_acc}}) {
                my $welded_id = $iworm_to_id{$welded_iworm_acc};
                print "\t$iworm_id" . "->" . "$welded_id [color=\"blue\"];\n";
            }
        }
    }
    
    if ($draw_type =~ /SCAFF|BOTH/) {
        ## add scaffolding links
        foreach my $iworm_acc (keys %scaffolds) {
            my $iworm_id = $iworm_to_id{$iworm_acc};
            unless (defined $iworm_id) { next; }
            
            foreach my $scaffolded_iworm_acc (keys %{$scaffolds{$iworm_acc}}) {
                my $scaffolded_id = $iworm_to_id{$scaffolded_iworm_acc};
                unless (defined $scaffolded_id) { next; }  # not all scaffolded entries were chosen for clustering if read support was insufficient
                
                print "\t$iworm_id" . "->" . "$scaffolded_id [color=\"green\"];\n";
            }
        }
    }
    
    print "}\n";

    exit(0);
    

}


####
sub parse_welds {
    my ($graph_file, $iworm_accs_href) = @_;
    
    my %welds;

    open (my $fh, $graph_file) or die "Error, cannot oepn file $graph_file";
    while (<$fh>) {
        chomp;
        if (/^\#Welding: >(a\d+;\d+)_\S+ to >(a\d+;\d+)_/) {
            my $iworm_acc_A = $1;
            my $iworm_acc_B = $2;

            ($iworm_acc_A, $iworm_acc_B) = sort ($iworm_acc_A, $iworm_acc_B);

            if ($iworm_accs_href->{$iworm_acc_A} || $iworm_accs_href->{$iworm_acc_B}) { # really should only have situations where they're both included
                
                $welds{$iworm_acc_A}->{$iworm_acc_B} = 1;
            }
        }
    }

    close $fh;

    return(%welds);
}


####
sub get_iworm_list_from_component {
    my ($graph_file, $component_no) = @_;

    my %iworm_accs;
    
    open (my $fh, $graph_file) or die "Error, cannot open file $graph_file";
    while (<$fh>) {
        if (/^>Component_$component_no /) {
            /\[iworm>(a\d+;\d+)_/ or die "Error, cannot parse iworm acc from $_ ";
            my $iworm_acc = $1;
            
            $iworm_accs{$iworm_acc} = 1;
        }
    }
    close $fh;


    return(%iworm_accs);
}


####
sub parse_scaffolded_pairs {
    my ($scaffolding_file, $iworm_accs_href) = @_;

    my %scaffolded_pairs;

    open (my $fh, $scaffolding_file) or die "Error, cannot open file $scaffolding_file";
    while (<$fh>) {
        chomp;
        my @x = split(/\t/);
        my ($iworm_acc_A, $iworm_acc_B) = ($x[0], $x[2]);
        
        ($iworm_acc_A, $iworm_acc_B) = sort ($iworm_acc_A, $iworm_acc_B);

        if ($iworm_accs_href->{$iworm_acc_A} || $iworm_accs_href->{$iworm_acc_B}) {
            
            $scaffolded_pairs{$iworm_acc_A}->{$iworm_acc_B} = 1;
        }
    }
    close $fh;

    return(%scaffolded_pairs);
}

