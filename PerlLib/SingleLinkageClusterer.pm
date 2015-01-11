package SingleLinkageClusterer;

## package not to be instantiated.  Just provides a namespace.

## Input: Array containing array-refs of pairs:
##               @_ = ( [1,2], [2,3], [6,7], [7,8], ...)
## Output: Array of all clusters as array-refs.
##              return ([1,2,3] , [6,7,8], ...)


use strict;
use warnings;


sub build_clusters {
    my @pairs = @_;
    my $pairfile = "$$.pairs";
    
    #must do mapping because cluster program doesn't like word chars, just ints.
    my %map_id_to_feat; 
    my %map_feat_to_id;
    my $id = 1;
    
    open (PAIRLIST, ">$pairfile") or die "Can't write $pairfile to /tmp";
    foreach my $pair (@pairs) {
        my ($a, $b) = @$pair;
        unless ($map_feat_to_id{$a}) {
            $map_feat_to_id{$a} = $id;
            $map_id_to_feat{$id} = $a;
            $id++;
        }
        unless ($map_feat_to_id{$b}) {
            $map_feat_to_id{$b} = $id;
            $map_id_to_feat{$id} = $b;
            $id++;
        }
        
        print PAIRLIST "$map_feat_to_id{$a} $map_feat_to_id{$b}\n";
    }
    close PAIRLIST;
    
    my $clusterfile = "$$.clusters";
    
    my $cluster_prog = "slclust"; 
        
    system "touch $clusterfile";
    unless (-w $clusterfile) { die "Can't write $clusterfile";}
    my $cmd = "$cluster_prog < $pairfile > $clusterfile";
    
    my $ret = system ($cmd);
    if ($ret) {
        die "ERROR: Couldn't run cluster properly via path: $cluster_prog.\ncmd: $cmd";
    }
    
    my @clusters;
    open (CLUSTERS, $clusterfile);
    
    while (my $line = <CLUSTERS>) {
        my @elements;
        while ($line =~ /(\d+)\s?/g) {
            push (@elements, $map_id_to_feat{$1});
        }
        if (@elements) {
            push (@clusters, [@elements]);
        }
    }
    
    close CLUSTERS;
    
    ## clean up
    unlink ($pairfile, $clusterfile);
    
    return (@clusters);
}


1;
