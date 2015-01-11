#!/usr/bin/env perl

use strict;
use warnings;


my $usage = "usage: $0 summary.dat [outprefix=argv[0]]\n\n";

my $file = $ARGV[0] or die $usage;
my $out_prefix = $ARGV[1] || $file;

main: {

    my %up_to_down;
    my %feature_to_sample_expr_val;

    open (my $fh, $file) or die $!;
    while (<$fh>) {
        chomp;
        my ($feature, $sampleA, $sampleB, $exprA, $exprB, $log_fold_change, $post_prob) = split(/\t/);

        
        $feature_to_sample_expr_val{$feature}->{$sampleA} = $exprA;
        $feature_to_sample_expr_val{$feature}->{$sampleB} = $exprB;
        
        my ($up, $down) = ($log_fold_change > 0) ? ($sampleA, $sampleB) : ($sampleB, $sampleA);
        
        $log_fold_change = abs($log_fold_change);
        
        $up_to_down{$feature}->{$up}->{$down} = { logFC => $log_fold_change,
                                                  Pval => $post_prob,
                                              };
        
        

    }
    

    my $class_up_outfile = "$out_prefix.class_up_priority";
    open (my $class_up_ofh, ">$class_up_outfile") or die "Error, cannot write to $class_up_outfile";
    
    my $class_down_outfile = "$out_prefix.class_down_priority";
    open (my $class_down_ofh, ">$class_down_outfile") or die "Error, cannot write to $class_down_outfile";
    
    my $graph_outfile = "$out_prefix.graph";
    open (my $graph_ofh, ">$graph_outfile") or die "Error, cannot write to $graph_outfile";
    
        
    my %up_cat_to_dat;
    my %up_class_counter;
    
    foreach my $feature (keys %up_to_down) {

        my ($top_up_list, $top_down_list) = &get_top_up_down_list($up_to_down{$feature}, 'up_priority');
        
        my $up_class = join(",", sort @$top_up_list);
        my $down_class = join(",", sort @$top_down_list);
        
        print $class_up_ofh join("\t", $feature, $up_class, $down_class) . "\n";
        


        &write_graph_entry($feature, $up_to_down{$feature}, $graph_ofh);
        
        my $up_class_expr = 0;
        
        my $min_up_expr = undef;
                
        foreach my $up ( @$top_up_list) {
            my $expr = $feature_to_sample_expr_val{$feature}->{$up};
            $up_class_expr += $expr;
            
            if ( (! defined $min_up_expr) || $min_up_expr > $expr) {
                $min_up_expr = $expr;
            }
        }
        $up_class_expr /= scalar(@$top_up_list);

        
        my $max_down_expr = 0;
        foreach my $down (@$top_down_list) {
            my $expr = $feature_to_sample_expr_val{$feature}->{$down};
            if ($expr > $max_down_expr) {
                $max_down_expr = $expr;
            }
        }
        

        my $delta_expr = $min_up_expr - $max_down_expr;
        my $priority = $delta_expr / ($max_down_expr + 1); # the 1 is a pseudocount to avoid ultra-small value comparisons.

        push (@{$up_cat_to_dat{$up_class}}, { up_expr => $up_class_expr,
                                              feature => $feature,
                                              up_class => $up_class,
                                              down_class => $down_class, 
                                              
                                              min_up_expr => $min_up_expr,
                                              max_down_expr => $max_down_expr,

                                              delta_expr => $delta_expr,
                                              
                                              priority =>  $priority,

                                          });
        
        

        $up_class_counter{$up_class}++;



        {
            ## similarly look for transcripts that are downregulated as hallmark features
            
            my ($top_up_list, $top_down_list) = &get_top_up_down_list($up_to_down{$feature}, 'down_priority');
            
            my $up_class = join(",", sort @$top_up_list);
            my $down_class = join(",", sort @$top_down_list);
            
            
            print $class_down_ofh join("\t", $feature, $up_class, $down_class) . "\n";
            
        }
        
    }
    
    close $class_up_ofh;
    close $class_down_ofh;
    close $graph_ofh;
    

    
    ## make prioritized list
    open (my $ofh_prioritized, ">$out_prefix.class_up_priority.ordered_by_expression") or die $!;

    foreach my $up_class (reverse sort {$up_class_counter{$a}<=>$up_class_counter{$b}} keys %up_cat_to_dat) {
        

        
        my @feature_structs = @{$up_cat_to_dat{$up_class}};

        @feature_structs = reverse sort {$a->{priority}<=>$b->{priority}} @feature_structs;
        

        my $num_features = scalar(@feature_structs);
        print $ofh_prioritized "## $up_class ($num_features)\n";

        foreach my $feature_struct (@feature_structs) {

            my $feature = $feature_struct->{feature};
            
            my @up_classes = split(/,/, $feature_struct->{up_class});
            my @down_classes = split(/,/, $feature_struct->{down_class});

            my @up_class_text;

            @up_classes = reverse sort {$feature_to_sample_expr_val{$feature}->{$a} <=> $feature_to_sample_expr_val{$feature}->{$b}} @up_classes;

            foreach my $up_class (@up_classes) {
                
                my $expr = sprintf("%.2f", $feature_to_sample_expr_val{$feature}->{$up_class});
                push (@up_class_text, "$up_class\($expr)");
            }
            
            my @down_class_text;
            @down_classes = reverse sort {$feature_to_sample_expr_val{$feature}->{$a} <=> $feature_to_sample_expr_val{$feature}->{$b}} @down_classes;
            foreach my $down_class (@down_classes) {
                my $expr = sprintf("%.2f", $feature_to_sample_expr_val{$feature}->{$down_class});
                push (@down_class_text, "$down_class\($expr)");
            }
            
            print $ofh_prioritized join("\t", $feature, join(",", @up_class_text), join(",", @down_class_text) . "\n");
        }
    }
    
    close $ofh_prioritized;
    

    exit(0);

}

####
sub get_top_up_down_list {
    my ($up_down_href, $priority_direction) = @_;

    my @structs;

    
    if ($priority_direction eq 'down_priority') {
        
        $up_down_href = &reverse_updown_list($up_down_href);
    }
    

    foreach my $up (keys %$up_down_href) {

        my $down_href = $up_down_href->{$up};
        
        my @down = keys %$down_href;

        my $struct = { up => $up,
                       down => [@down],
                       num => scalar @down,
                   };

        push (@structs, $struct);
    }

    @structs = reverse sort {$a->{num}<=>$b->{num}} @structs;

    my $top_struct = shift @structs;
    my @top_structs = ($top_struct);
    
    my $top_num = $top_struct->{num};
    
    while (@structs) {
        my $struct = shift @structs;
        if ($struct->{num} == $top_num) {
            push (@top_structs, $struct);
        }
        else {
            last;
        }
    }
    
    my @top;
    my %bottom;
    foreach my $struct (@top_structs) {
        push (@top, $struct->{up});
        
        foreach my $down (@{$struct->{down}}) {
            $bottom{$down}++;
        }
    }
    
    @top = sort @top;
    my @bottom = sort keys %bottom;

    if ($priority_direction eq 'down_priority') {
        
        # switch them around.
        
        my @orig_bottom = @bottom;
        @bottom = @top;
        @top = @orig_bottom;
        
    }


    return(\@top, \@bottom);
}


####
sub write_graph_entry {
    my ($feature, $graph_href, $graph_ofh) = @_;
    
    my @nodes;
    
    foreach my $up_sample (keys %$graph_href) {

        foreach my $down_sample (keys %{$graph_href->{$up_sample}}) {

            my $struct = $graph_href->{$up_sample}->{$down_sample};
            
            my $logFC = $struct->{logFC};
            
            push (@nodes, "$up_sample,$down_sample,$logFC");
        }

    }
    print $graph_ofh join("\t", $feature, @nodes) . "\n";
    
    return;
}
    


####
sub reverse_updown_list {
    my ($up_down_href) = @_;

    my %down_up_href;

    foreach my $up (keys %$up_down_href) {

        foreach my $down (keys %{$up_down_href->{$up}}) {

            $down_up_href{$down}->{$up} = $up_down_href->{$up}->{$down};
        }
    }


    return(\%down_up_href);

}

