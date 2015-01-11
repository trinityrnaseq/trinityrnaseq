#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;

my $usage = "usage: $0 file.paralog_clusters file.maps\n\n";

my $paralogs_file = $ARGV[0] or die $usage;
my $maps_file = $ARGV[1] or die $usage;


main: {

	my @paralog_clusters = &get_paralog_clusters($paralogs_file);
	
	my %FL;
	open (my $fh, $maps_file) or die "Error, cannot open file $maps_file";
	while (<$fh>) {
		chomp;
		my @x = split(/\t/);
		my $genes_list = $x[0];
		foreach my $gene (split(/,/, $genes_list) ) {
			$FL{$gene} = 1;
            my ($trans_id, $gene_id) = split(/;/, $gene);
            if ($gene_id) {
                $FL{$gene_id} = 1;
            }
        }
	}
	close $fh;
	

	## stats to capture
	my $number_FL_paralogs = 0;
	my $number_FL_clusters_in_entirety = 0;
	my $number_paralog_clusters_any = 0;

	my $total_paralogs = 0;
	my $total_clusters = scalar(@paralog_clusters);
	my $total_clusters_at_least_2 = 0;

	foreach my $paralog_cluster (@paralog_clusters) {
		
		my $missing = 0;
		my $got_one = 0;
		foreach my $paralog (@$paralog_cluster) {
			$total_paralogs++;
			if ($FL{$paralog}) {
				$number_FL_paralogs++;
				$got_one++;
				print "$paralog\[YES]\t";
			}
			else {
				$missing = 1;
				print "$paralog\[NO]\t";
			}
		}
		print "\n";
		if (! $missing) {
			$number_FL_clusters_in_entirety++;
		}
		if ($got_one) {
			$number_paralog_clusters_any++;
			if ($got_one > 1) {
				$total_clusters_at_least_2++;
			}
		}
	}
		
	print "\n\n";
	print "$number_FL_paralogs FL paralogs / $total_paralogs = " . sprintf("%.2f", $number_FL_paralogs/$total_paralogs*100) . "\n";
	print "$number_paralog_clusters_any clusters with at least one FL entry / $total_clusters = " . sprintf("%.2f", $number_paralog_clusters_any / $total_clusters * 100) . "\n";
	
	print "$total_clusters_at_least_2 clusters with at least 2 FL entries / $total_clusters = " . sprintf("%.2f", $total_clusters_at_least_2 / $total_clusters * 100) . "\n";
	print "$number_FL_clusters_in_entirety clusters fully covered by FL entries / $total_clusters = " . sprintf("%.2f", $number_FL_clusters_in_entirety / $total_paralogs * 100) . "\n";
	

}


####
sub get_paralog_clusters {
    my ($paralogs_file) = @_;
    
	my @paralog_clusters;
	
	open (my $fh, $paralogs_file) or die "Error, cannot find paralog info";
	while (<$fh>) {
		chomp;
		my @eles = split(/\s+/);
		push (@paralog_clusters, [@eles]);
	}
	close $fh;

	return(@paralog_clusters);
}

