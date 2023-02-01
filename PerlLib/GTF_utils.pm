package GTF_utils;

use strict;
use warnings;
use Gene_obj;
use Gene_obj_indexer;
use GTF;
use Carp;
use Data::Dumper;
use Overlap_piler;

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

    my %seen;
    
    for my $gene_obj (@$gene_objs) {

        my $gene_id = $gene_obj->{TU_feat_name};
        

        if ($seen{$gene_id}) {
            confess "Error, already processed gene: $gene_id\n"
                . " here: " . $gene_obj->toString() . "\n"
                . " and earlier: " . $seen{$gene_id}->toString();
            
        }
        
        $seen{$gene_id} = $gene_obj;


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
    my %noncoding_features;

    my %gene_id_to_source;
    my %gene_id_to_name;

    my %gene_id_to_seq_name;
    my %gene_id_to_gene_name;

    my %gene_id_to_gene_type;
    my %transcript_id_to_transcript_type;

    
    my %coding_genes;
    
        
    open (my $fh, $gtf_filename) or die "Error, cannot open $gtf_filename";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } # comment line.
        
        chomp;
        my ($seqname, $source, $type, $lend, $rend, $score, $strand, $gtf_phase, $annot) = split (/\t/);
        
        my ($end5, $end3) = ($strand eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $annot =~ /gene_id \"([^\"]+)\"/  or confess "Error, cannot get gene_id from $annot of line\n$_";
        my $gene_id = $1;

        if (my $sn = $gene_id_to_seq_name{$gene_id}) {
            if ($sn ne $seqname) {
                $gene_id = "$seqname" . "|" . "$gene_id"; # make unique per scaffold.
            }
        }
        $gene_id_to_seq_name{$gene_id} = $seqname;
        
        $gene_id_to_source{$gene_id} = $source;


        if ($annot =~ /name \"([^\"]+)\"/) {
            my $name = $1;
            $gene_id_to_name{$gene_id} = $name;
        }

        my $gene_name = "";
        if ($annot =~ /gene_name \"([^\"]+)\"/) {
            $gene_name = $1;
            $gene_id_to_gene_name{$gene_id} = $gene_name;
        }

        if ($annot =~ /gene_type \"([^\"]+)\"/) {
            my $gene_type = $1;
            $gene_id_to_gene_type{$gene_id} = $gene_type;
        }

        
		# print "gene_id: $gene_id, transcrpt_id: $transcript_id, $type\n";

        if ($type eq 'transcript' || $type eq 'gene') { next; } # capture by exon coordinates

        my $transcript_id;
        if ($annot =~ /transcript_id \"([^\"]+)\"/) {
            $transcript_id = $1;
            
            if ($annot =~ /transcript_type \"([^\"]+)\"/) {
                my $transcript_type = $1;
                $transcript_id_to_transcript_type{$transcript_id} = $transcript_type;
            }
        }
        else {
            print STDERR "Skipping line: $_, no transcript_id value provided\n";
            next;
        }
        
        if ($type eq 'CDS' || $type eq 'stop_codon' || $type eq 'start_codon') {
            push (@{$gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{CDS}}, [$end5, $end3] );
            push (@{$gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}}, [$end5, $end3] );

            $coding_genes{$gene_id}++;
        }
        elsif ($type eq "exon" || $type =~ /UTR/) {
            push (@{$gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}}, [$end5, $end3] );
        }
        elsif ($type =~ /Selenocysteine/) {
            # no op
        }
        else {
            ## assuming noncoding feature
            push (@{$noncoding_features{$seqname}->{$type}->{$gene_id}->{$transcript_id}}, [$end5, $end3] );
        }

    }
    close $fh;
    

    ## create gene objects.
 
    my @top_gene_objs;

    my %seen;
    foreach my $seqname (keys %gene_transcript_data) {


        {
            ##################################
            ## Process protein-coding genes:
            
            my $genes_href = $gene_transcript_data{$seqname};
            
            foreach my $gene_id (keys %$genes_href) {
             
                if ($seen{$gene_id}) {
                    print STDERR ("Error, already saw $gene_id,$seqname as $gene_id,$seen{$gene_id}\nSkipping.");
                    next;
                }
                $seen{$gene_id} = $seqname;
                
                my $transcripts_href = $genes_href->{$gene_id};
                
                my $source = $gene_id_to_source{$gene_id};
                
                my @gene_objs;
                
                foreach my $transcript_id (keys %$transcripts_href) {
                    
                    my $coord_types_href = $transcripts_href->{$transcript_id};
                    
                    my $CDS_coords_aref = $coord_types_href->{CDS};
                    my $mRNA_coords_aref = $coord_types_href->{mRNA};
                    
                    
                    #print STDERR "Before, CDS: " . Dumper($CDS_coords_aref);
                    #print STDERR "Before, exons: " . Dumper($mRNA_coords_aref);
                    
                    
                    my $CDS_coords_href = &_join_overlapping_coords($CDS_coords_aref);
                    my $mRNA_coords_href = &_join_overlapping_coords($mRNA_coords_aref);
                    
                    #print STDERR "CDS: " . Dumper($CDS_coords_href);
                    #print STDERR "mRNA: " . Dumper ($mRNA_coords_href);
                    
                    
                    my $gene_obj = new Gene_obj();
                    $gene_obj->populate_gene_object($CDS_coords_href, $mRNA_coords_href);
                    
                    $gene_obj->{TU_feat_name} = $gene_id;
                    $gene_obj->{Model_feat_name} = $transcript_id;
                    if (my $name = $gene_id_to_name{$gene_id}) {
                        $gene_obj->{com_name} = $name;
                    }
                    else {
                        $gene_obj->{com_name} = $transcript_id;
                    }
                    if (my $gene_name = $gene_id_to_gene_name{$gene_id}) {
                        $gene_obj->{gene_name} = $gene_name;
                    }
                    $gene_obj->{asmbl_id} = $seqname;
                    $gene_obj->{source} = $source;
                    
                    if (my $gene_type = $gene_id_to_gene_type{$gene_id}) {
                        $gene_obj->{gene_type} = $gene_type;
                    }
                    if (my $transcript_type = $transcript_id_to_transcript_type{$transcript_id}) {
                        $gene_obj->{transcript_type} = $transcript_type;
                    }
                    
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
        
        
        {
            ################################
            ## Process noncoding features ##
            ################################
            
            my $ncgene_types_href = $noncoding_features{$seqname};
            
            if (ref $ncgene_types_href) {
                
                foreach my $nc_type (keys %$ncgene_types_href) {
                    
                    my $gene_ids_href = $ncgene_types_href->{$nc_type};
                    foreach my $gene_id (keys %$gene_ids_href) {

                        if (exists $coding_genes{$gene_id}) {
                            print STDERR "Warning: Skipping $gene_id ($nc_type) as this gene is already included as a coding gene.\n";
                            next;
                        }

                        my $trans_ids_href = $gene_ids_href->{$gene_id};
                        foreach my $trans_id (keys %$trans_ids_href) {
                            
                            my @coordsets = @{$trans_ids_href->{$trans_id}};
                            my %coords;
                            foreach my $coordset (@coordsets) {
                                my ($end5, $end3) = @$coordset;
                                $coords{$end5} = $end3;
                            }

                            my $gene_obj = new Gene_obj();
                            
                            $gene_obj->populate_gene_object({}, \%coords);
                            $gene_obj->{asmbl_id} = $seqname;
                            $gene_obj->{TU_feat_name} = $gene_id;
                            $gene_obj->{Model_feat_name} = $trans_id;
                            $gene_obj->{gene_type} = $nc_type;
                            if (my $name = $gene_id_to_name{$gene_id}) {
                                $gene_obj->{com_name} = $name;
                            }
                            else {
                                $gene_obj->{com_name} = $trans_id;
                            }
                            
                            
                            #print STDERR $gene_obj->toString();
                           

                            push (@top_gene_objs, $gene_obj);
                        }
                    }
                }
            }
        }
    }
    return (\@top_gene_objs);
}



####
sub _join_overlapping_coords {
    my $coords_aref = shift;

    unless (ref $coords_aref) {
        return ({});
    }

    my $orient;
    
    my @coords;
    
    my $inferred_orient;

    foreach my $coordset (@$coords_aref) {
        my ($end5, $end3) = @$coordset;
        
        my $orient;
        if ($end5 < $end3) {
            $orient = '+';
        }
        elsif ($end5 > $end3) {
            $orient = '-';
        
        }
        
        if ( (! defined $inferred_orient) && defined $orient) {
            $inferred_orient = $orient;
        }
        elsif ( defined($orient) && $orient ne $inferred_orient) {
            die "Error, conflicting orientation info: " . Dumper($coords_aref);
        }
        
        my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
        
        push (@coords, [$lend, $rend]);
    }
    
    #print STDERR "coords: " . Dumper(@coords);

    my @piles = &Overlap_piler::simple_coordsets_collapser(@coords);

    #print STDERR "piles: " . Dumper(@piles);

    @piles = sort {$a->[0] <=> $b->[0]} @piles;
    
    ## join adjacent piles
    my @joined_piles = shift @piles;
    foreach my $pile (@piles) {
        my ($pile_lend, $pile_rend) = @$pile;
        if ($pile_lend == $joined_piles[$#joined_piles]->[1] + 1) {
            # adjacent
            $joined_piles[$#joined_piles]->[1] = $pile_rend;
        }
        else {
            push (@joined_piles, $pile);
        }
    }

    
    my %new_coords;
    foreach my $pile (@joined_piles) {
        my ($lend, $rend) = @$pile;
        my ($end5, $end3) = ($inferred_orient eq '+') ? ($lend, $rend) : ($rend, $lend);
        $new_coords{$end5} = $end3;
    }


    return(\%new_coords);
}
    
    
            
1; #EOM
