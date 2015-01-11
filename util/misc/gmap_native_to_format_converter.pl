#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../PerlLib");
use Gene_obj;

my $usage = "usage: $0 file.gmap (BED|GTF)\n\n";


my $gmap_file = $ARGV[0] or die $usage;
my $out_format = $ARGV[1] or die $usage;

unless ($out_format =~ /^(BED|GTF)$/) {
    die $usage;
}


my %data;

open (my $fh, $gmap_file) or die $!;
while (<$fh>) {
    if (/^>(\S+)/) {
        if (%data) {
			&format_output(%data);
		}
		
		%data = ();
		
		my $transcript_acc = $1;
    
		%data = (acc => $transcript_acc,
				 segments => {},  # end5 => end3
			);
		
	}
    elsif (/\s*([\+\-])([^\:]+):(\d+)-(\d+)\s+\((\d+)-(\d+)\)\s+(\d[^\%]+\%)/) {
        #print " $transcript_acc\t$_";
        my $orient = $1;
        my $genome_acc = $2;
        my $genome_end5 = $3;
        my $genome_end3 = $4;
        my $transcript_end5 = $5;
        my $transcript_end3 = $6;
        my $percent_identity = $7;
        
		
		$data{segments}->{$genome_end5} = $genome_end3;
		$data{scaff} = $genome_acc;
	}
}

if (%data) {
	&format_output(%data);
}

exit(0);


####
sub format_output {
	my %data = @_;

	my $gene_obj = new Gene_obj();
	
	my $coords_href = $data{segments};
	my $acc = $data{acc};
	my $genome_scaff = $data{scaff};

	if (%$coords_href) {
		
		$gene_obj->populate_gene_object($coords_href, $coords_href);
		
		$gene_obj->{com_name} = $acc;
		$gene_obj->{asmbl_id} = $genome_scaff;
		$gene_obj->{TU_feat_name} = "g|$acc";
        $gene_obj->{Model_feat_name} = "t|$acc";
        
        if ($out_format eq "BED") {
            print $gene_obj->to_BED_format();
        }
        elsif ($out_format eq "GTF") {
            print $gene_obj->to_transcript_GTF_format() . "\n";
        }
        else {
            die "Error, cannot process format $out_format";
            # shouldn't get here if proper option values processed at top.
        }
        

	}
	
	return;
}
