#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use Gene_obj;
use Fasta_reader;

my $usage = "usage: $0 transcripts.gtf genome.fasta\n\n";

my $gtf_file = $ARGV[0] or die $usage;
my $genome_fasta_file = $ARGV[1] or die $usage;

main: {
	

    my $fasta_reader = new Fasta_reader($genome_fasta_file);
    my %genome_seqs = $fasta_reader->retrieve_all_seqs_hash();
    

	my %genome_trans_to_coords;
	
	open (my $fh, $gtf_file) or die "Error, cannot open file $gtf_file";
	while (<$fh>) {
		chomp;
		
		unless (/\w/) { next; }
		
		my @x = split(/\t/);
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];

		my $orient = $x[6];
		
		my $info = $x[8];
		
		unless ($type eq 'exon') { next; }

		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/\s+/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{gene_id} or die "Error, no gene_id at $_";
		my $trans_id = $atts{transcript_id} or die "Error, no trans_id at $_";
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$genome_trans_to_coords{$scaff}->{$gene_id}->{$trans_id}->{$end5} = $end3;

	}


	## Output genes in gff3 format:

    

	foreach my $scaff (sort keys %genome_trans_to_coords) {

		my $genes_href = $genome_trans_to_coords{$scaff};

        my $genome_seq = $genome_seqs{$scaff} or die "Error, cannot find genome sequence for acc: $scaff";
        

		foreach my $gene_id (keys %$genes_href) {

			my $trans_href = $genes_href->{$gene_id};

            my $intron_text = "";
            
			foreach my $trans_id (keys %$trans_href) {

				my $coords_href = $trans_href->{$trans_id};

				my $gene_obj = new Gene_obj();

				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = $trans_id;
				$gene_obj->{com_name} = "$gene_id $trans_id";
				
				$gene_obj->{asmbl_id} = $scaff;
				
				$gene_obj->populate_gene_object($coords_href, $coords_href);
			
                my $orientation = $gene_obj->get_orientation();

				#print $gene_obj->to_BED_format();
                
                my @intron_coords = $gene_obj->get_intron_coordinates();
                
                foreach my $intron (@intron_coords) {
                    my ($end5, $end3) = @$intron;
                    my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
                    my $left_splice_dinuc = substr($genome_seq, $lend-1, 2);
                    my $right_splice_dinuc = substr($genome_seq, $rend-1-1, 2);
                    
                    $intron_text .= join("\t", $gene_id, $trans_id, $scaff, "$lend-$rend", $orientation, "$left_splice_dinuc..$right_splice_dinuc") . "\n";
                }
            }
                
            if ($intron_text) {
                print "$intron_text\n";
            }
			            
		}
	}


	exit(0);
}

