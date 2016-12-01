#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::RealBin/../../../PerlLib");

use SingleLinkageClusterer;
use PSL_parser;
use Process_cmd;

use Getopt::Long qw(:config no_ignore_case bundling);

use DelimParser;
use File::Basename;


my $usage = <<__EOUSAGE__;

####################################################################################
#
#  --target <string>
#  --query <string>
#
#  --min_per_id <int>          default(95)
#  --max_per_gap <int>    default(5)
#  --min_per_length <float>   default(100)
#
#  --forward_orient 
#
#  --allow_non_unique_mappings    a single query transcript counts towards all 
#                                   transcripts encapsulated
#
#  --reuse                       don't prompt for reusing blat if the output already exists.
#  --no_reuse                    dont prompt and do NOT reuse existing outputs
#
#  --out_prefix <string>         output prefix
#
#####################################################################################

__EOUSAGE__

    ;

my $target;
my $query;

my $min_per_id = 95;
my $max_per_gap = 5;
my $min_per_length = 100;

my $forward_orient = 0;

my $help_flag;

my $allow_non_unique_mappings = 0;

my $reuse_flag = 0;
my $no_reuse_flag = 0;

my $out_prefix = "";

&GetOptions ( 'h' => \$help_flag,
              'target=s' => \$target,
              'query=s' => \$query,
              'min_per_id=i' => \$min_per_id,
              'forward_orient' => \$forward_orient,
              'max_per_gap=i' => \$max_per_gap,
              'allow_non_unique_mappings' => \$allow_non_unique_mappings,
              'min_per_length=f' => \$min_per_length,
	      
              'reuse' => \$reuse_flag,
              'no_reuse' => \$no_reuse_flag,

              'out_prefix=s' => \$out_prefix,

              );

unless ($target && $query) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}



my %mapping_info;


main: {

    unless ($out_prefix) {
        $out_prefix = basename($query);
    }
    
	my $blat_output = "$out_prefix.pslx";
	
    my  $blat_prog = `which blat`;
    chomp $blat_prog;
    unless ($blat_prog) {
        die "Error, cannot find a blat program in your path.\n";
    }
    
	my $cmd = "$blat_prog -t=dna -q=dna -out=pslx $target $query $blat_output > /dev/null";
	
	my $run_blat = 0;
	
	if (-s $blat_output) {
        
        if ($no_reuse_flag) {
            $run_blat = 1;
        }
        elsif (! $reuse_flag) {
            print STDERR "Regenerate blat output file?";
            my $response = <STDIN>;
            if ($response =~ /^y/i) {
                
                $run_blat = 1;
            }
        }
    }
	else {
		$run_blat = 1;
	}
    
	if ($run_blat) {
		
		my $ret = system($cmd);
		if ($ret) {
			die "Error, cmd: $cmd died with ret $ret";
		}
		else {
			print STDERR "Blat done.  Parsing output...\n";
		}
	}
	
	my %gene_to_encapsulating_assemblies = &examine_blat_mappings($blat_output);
	
	my @pairs;
	foreach my $gene (keys %gene_to_encapsulating_assemblies) {
		
		my @assemblies = @{$gene_to_encapsulating_assemblies{$gene}};
        
        if ($allow_non_unique_mappings) {
            $gene =~ s/GENE://;
            foreach my $asm (@assemblies) {
                $asm =~ s/ASM://;
            }
            
            print "$gene\t" . join(",", @assemblies) . "\n";
            next;
        }
		

		foreach my $assembly (@assemblies) {
			
			push (@pairs, [$gene, $assembly]);
			
		}
	}

    if ($allow_non_unique_mappings) {
        ## done.
        exit(0);
    }
    

    ## Find best transcript mapping to gene for those that are full-length
    
    #print "Input as pairs: @pairs\n";

	my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);

    #print "Clusters: @clusters\n";
    #die;

	my $selected_blat_file = "$blat_output.FL_selected";
	open (my $ofh, ">$selected_blat_file") or die $!;
	

	foreach my $cluster (@clusters) {
		my @genes = grep { /GENE:/ } @$cluster;
		my @asms = grep { /ASM:/ } @$cluster;
		
		foreach my $gene (@genes) {
			$gene =~ s/GENE://;
		}
		foreach my $asm (@asms) {
			$asm =~ s/ASM://;
		}
		
        ## find best matching, resolve overlaps
			
        &resolve_overlaps(\@genes, \@asms, $ofh);
		

	}
	
	exit(0);


}

####
sub resolve_overlaps {
	my ($genes_aref, $asms_aref, $ofh) = @_;

	my @genes = @$genes_aref;
	my @asms = @$asms_aref;
	
	my @mappings;

	foreach my $asm (@asms) {
		
		foreach my $gene (@genes) {
			
			my $key = join("$;", $gene, $asm);
			if (exists $mapping_info{$key}) {
				
				my $struct = $mapping_info{$key};
				push (@mappings, $struct);
			}
		}
	}

	@mappings = reverse sort {$a->{score}<=>$b->{score}} @mappings;
	
	
	my @report;
	
	foreach my $mapping (@mappings) {
		# ensure no overlap among assembly and gene coordinates
		if (! &overlaps_asm($mapping, \@report)) {
			push (@report, $mapping);
		}
	}
	
	
	## regroup annotations and assemblies

	my %asm_links;
	foreach my $entry (@report) {
		my ($gene_id, $trans_id) = ($entry->{gene_id}, $entry->{trans_id});
		push (@{$asm_links{$trans_id}}, $gene_id);
		
		my $key = join("$;", $gene_id, $trans_id);
		my $struct = $mapping_info{$key} or die "Error, cannot find alignment for $key";
		
		print $ofh $struct->{line} . "\n";
        
	}
	
	foreach my $asm (keys %asm_links) {
		my @genes = @{$asm_links{$asm}};
        
		print join(",", @genes) . "\t$asm\n";
	}
    
	return;
}


####
sub overlaps_asm {
	my ($mapping, $report_aref) = @_;

	my ($map_lend, $map_rend) = ($mapping->{trans_lend}, $mapping->{trans_rend});
	my $gene_id = $mapping->{gene_id};
	my $trans_id = $mapping->{trans_id};

	foreach my $entry (@$report_aref) {

		my $entry_trans_id = $entry->{trans_id};
		my $entry_gene_id = $entry->{gene_id};

		if ($gene_id eq $entry_gene_id) {
			# already used gene
			return(1);
		}

		elsif ($entry_trans_id eq $trans_id) {
		
			my ($lend, $rend) = ($entry->{trans_lend}, $entry->{trans_rend});
			
			
			if ($lend < $map_rend && $rend > $map_lend) {
				return(1);
			}
		}
	}
	
	return(0); # no meaningful overlap
}



####
sub examine_blat_mappings {
	my ($blat_output) = @_;

	my %mappings;

    my $cmd = "$FindBin::Bin/blat_psl_to_align_summary_stats.pl $blat_output > $blat_output.stats";
    &process_cmd($cmd);

    open(my $fh, "$blat_output.stats") or die "Error, cannot open file $blat_output.stats";
    my $tab_reader = new DelimParser::Reader($fh, "\t");

    while (my $row = $tab_reader->get_row()) {
		
        #print $psl_entry->toString();

        my $line = $tab_reader->reconstruct_line_from_row($row);
        
		my $trans_assembly = $row->{query};
		my $gene_id = $row->{target};
		my $per_id = $row->{per_id};
		
        my $orient = $row->{strand};
        
		if ($forward_orient && $orient eq '-') { next; }
		if ($per_id < $min_per_id) { next; }

		my $percent_gapped = $row->{per_gap};
        
		if ($percent_gapped > $max_per_gap ) {
			#print STDERR "\%gapped = $percent_gapped, so skipping.\n";
            next;  # too gappy
		}
		

		my ($gene_lend, $gene_rend) = ($row->{target_lend}, $row->{target_rend});

		my ($trans_lend, $trans_rend) = ($row->{query_lend}, $row->{query_rend});

		my $pct_len = $row->{per_len};
        
		#print STDERR "PERCENT_LEN: $pct_len (vs. $min_per_length)\n";

		if ($pct_len >= $min_per_length) {
		    
		    #print STDERR "OK!!!\n";
		    
			my $pair = join("$;", $gene_id, $trans_assembly);

			
			## score the alignment:
			my $score = $row->{score};

			
			my $struct = $mapping_info{$pair};
			
			if ( (! $struct) || $score > $struct->{score}) {
				
				push (@{$mappings{"GENE:$gene_id"}}, "ASM:$trans_assembly");
								
				$mapping_info{$pair} = { "per_id" => $per_id,
										 
										 
										 gene_id => $gene_id,
										 trans_id => $trans_assembly,
										 
										 gene_lend => $gene_lend,
										 gene_rend => $gene_rend,
										 
										 trans_lend => $trans_lend,
										 trans_rend => $trans_rend,
										 
										 length => $gene_rend - $gene_lend + 1,
										 per_id => $per_id,
										 score => $score,

										 strand => $orient,
                                         
                                         line => $line,
                                         
                };
			}
			            
		}
	}
    
	return(%mappings);
}







####
sub parse_sequence_lengths {
	my ($seqlen_file) = @_;

	my %acc_to_len;

	open (my $fh, $seqlen_file) or die "Error, cannot open file $seqlen_file";
	while (<$fh>) {
		chomp;
		my @x = split(/\t/);
		my $acc = $x[1];
		my $len = $x[0];
		
		$acc_to_len{$acc} = $len;
	}

	close $fh;
		
	return(%acc_to_len);
}
