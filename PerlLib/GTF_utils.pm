package GTF_utils;

use strict;
use warnings;
use Gene_obj;
use Gene_obj_indexer;
use GTF;
use Carp;
use Data::Dumper;

####
sub index_GTF_gene_objs {
  my ($gtf_filename, $gene_obj_indexer) = @_;  ## gene_obj_indexer can be a simple  hashref {}
  
  unless ($gtf_filename && $gene_obj_indexer) {
	  confess "Error, need gtf_filename and gene_obj_indexer as perams";
  }

  return index_GTF_gene_objs_from_GTF($gtf_filename,$gene_obj_indexer);
}

sub index_GTF_gene_objs_from_GTF {
    my ($gtf_filename, $gene_obj_indexer) = @_;
    
    ## 
    #print STDERR "\n-caching genes.\n";
    my %seqname_map;

    my $gene_objs = GTF_to_gene_objs($gtf_filename);
    
    for my $gene_obj (@$gene_objs) {

        my $gene_id = $gene_obj->{TU_feat_name};
        my $seqname = $gene_obj->{asmbl_id};
        
		if (ref $gene_obj_indexer eq "HASH") {
			$gene_obj_indexer->{$gene_id} = $gene_obj;
		}
		else {
			$gene_obj_indexer->store_gene($gene_id, $gene_obj) if(ref $gene_obj_indexer);
		}
		
        # add to gene list for asmbl_id
        my $gene_list = $seqname_map{$seqname};
        unless (ref $gene_list) {
            $gene_list = $seqname_map{$seqname} = [];
        }
        push (@$gene_list, $gene_id);
    }
    return (\%seqname_map);
}

sub GTF_to_gene_objs {
    my ($gtf_filename) = @_;

    my %gene_transcript_data;

    my %gene_id_to_source;

    open (my $fh, $gtf_filename) or die "Error, cannot open $gtf_filename";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } # comment line.
        
        chomp;
        my ($seqname, $source, $type, $lend, $rend, $score, $strand, $gtf_phase, $annot) = split (/\t/);
        
        my ($end5, $end3) = ($strand eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $annot =~ /gene_id \"([^\"]+)\"/  or confess "Error, cannot get gene_id from $annot of line\n$_";
        my $gene_id = $1;
        
        $gene_id_to_source{$gene_id} = $source;

        $annot =~ /transcript_id \"([^\"]+)\"/  or confess "Error, cannot get transcript_id from $annot of line\n$_";
        my $transcript_id = $1;

		# print "gene_id: $gene_id, transcrpt_id: $transcript_id, $type\n";
        
        if ($type eq "exon" || $type eq 'CDS' || $type eq 'stop_codon') {
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{CDS}->{$end5} = $end3;
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}->{$end5} = $end3;
        }
        if ($type =~ /UTR/) {
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}->{$end5} = $end3;
        }
    }
    close $fh;
    
    ## create gene objects.
 
    my @top_gene_objs;
    
    foreach my $seqname (keys %gene_transcript_data) {

        my $genes_href = $gene_transcript_data{$seqname};
        
        foreach my $gene_id (keys %$genes_href) {
            
            my $transcripts_href = $genes_href->{$gene_id};
            
            my $source = $gene_id_to_source{$gene_id};

            my @gene_objs;

            foreach my $transcript_id (keys %$transcripts_href) {

                my $coord_types_href = $transcripts_href->{$transcript_id};

                my $CDS_coords_href = $coord_types_href->{CDS};
                my $mRNA_coords_href = $coord_types_href->{mRNA};
                

                my $gene_obj = new Gene_obj();
                $gene_obj->populate_gene_object($CDS_coords_href, $mRNA_coords_href);
                
                $gene_obj->{TU_feat_name} = $gene_id;
                $gene_obj->{Model_feat_name} = $transcript_id;
                $gene_obj->{com_name} = $transcript_id;
                $gene_obj->{asmbl_id} = $seqname;
                $gene_obj->{source} = $source;
                
                $gene_obj->join_adjacent_exons();

                push (@gene_objs, $gene_obj);
            }
			
                    
            ## want single gene that includes all alt splice variants here
            if(scalar(@gene_objs)) {
                my $template_gene_obj = shift @gene_objs;
                foreach my $other_gene_obj (@gene_objs) {
                    $template_gene_obj->add_isoform($other_gene_obj);
                }
                push (@top_gene_objs, $template_gene_obj);       
            
                # print $template_gene_obj->toString(); 
                
                
            }
        }
    }
    return (\@top_gene_objs);
}


1; #EOM
