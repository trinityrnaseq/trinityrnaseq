#!/usr/bin/env perl

use warnings;

#CHANGELOG
#05Jun09 - 	Added option to ignore X/Ns when estimating size of Query 
#			(to allow blasting specific strains of an assembly only after convert_project and replacing @s with Xs)
#			Provide the initial fastafile used to blast or 
#11Jun09 -	Added option to extract hits and misses ('CDS' and 'UTR')
#
#
#
#
our $VERSION = '1.5';

####
###	TODO: add dbname on the top entry of the hashes so that we can have multiple dbs in same file?
# TODO: implement DBM to store hash http://www.perl.com/pub/a/2006/02/16/mldbm.html

=head1 NAME 

analyze_assembly_blast.pl - Provides some simple descriptive stats in a SearchIO BLAST report.
Useful for example to investigate the coverage of an experiment-derived FASTA file A versus a reference FASTA file B.

=head1 VERSION

 Version 1.4
 
=head1 USAGE

 analyze_assembly_blast.pl [-cs <bitscore cutoff> -ce <evalue cutoff> -l <limit top hits>] <one or more BLAST files> 

	-cs|--cutoff-score:i	Specify cutoff score for report. Defaults to 80
	-ce|--cutoff-evalue:s	Specify cutoff evalue for report. Defaults to 1e-5
	-f|--format:s		Specify format (currently only BLAST supported)
	-d			Debug output. -dd is more detailed (only small reports: high memory usage)
	-l|--limit:i		Limit processing to these many top hits for each query. Will also limit the hash output (so cannot reuse a hash for a larger limit)
	-uh|--use_hash		If for some reason you want to redo the search (e.g. change cutoff, reduce limit etc), 
				this switch will speed things up significantly. No need to provide the original BLAST file. Can be used for more stringent -cs/-ce than the one used to build the hash, but clearly, cannot use a less stringent -cs/-ce
	-nt|--notimer		Specify to switch off timer.
	--lines:i		Specify number of lines in hash (to speed up timer stats)
	-x|maskxn		Don't include Ns or Xs when calculating query sequence length. Give FASTA used to BLAST or a number to set size of query.
	-extract:s		If you give the FASTA used as input to BLAST, I can give you the regions matching and not matching as FASTA. Cannot be used with -uh

=head1 AUTHORS

	Alexie Papanicolaou 1 2
	with help from Paul Wilkinson 2

	1 Max Planck Institute for Chemical Ecology, Germany
	2 Centre for Ecology and Conservation, University of Exeter, UK
	alexie@butterflybase.org
	
=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. 
Please note that incorporating the whole software or parts of its code in proprietary software 
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

Unfortunately tiles hsps only works with normal blast report, not tabulated output...
also, tblastx is a bit tricky... not sure if tiled works properly there (blastn is fine)

Verification from Jason S and Chris Field. BioPerl tiling is crap, best to use wu-blast apparently
with the -links option.

=cut


use strict;
use warnings;
use Pod::Usage;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Search::SearchUtils;
use Bio::Index::Fasta;
use Statistics::Descriptive;
#use Time::Progress;
#my $timer = new Time::Progress;
use Getopt::Long;
$| = 1;

# Declare options
my (
	 $debug,    $debug2, 
	 $limit,     $use_hash, $multiple_blasts,
	 $debugfile, $notimer,  $no_check
);
my $store_hash   = 1;         # forced because queries depend on it.
my $cut_score    = 80;
my $cut_evalue   = '1e-5';
my $report_style = "blast";
my $hash_lines;

# other global variables
my ( $query_tlength, $query_tlength_with_hits,$idfile,%ids,$maskxn,$extract,$extract_inx,%extr_hash,$single_copy );
GetOptions(
	'cs|cutoff-score:i'  => \$cut_score,
	'ce|cutoff-evalue:s' => \$cut_evalue,
	'f|format:s'         => \$report_style,
	'd'                  => \$debug,
	'dd'                 => \$debug2,
	'l|limit:i'          => \$limit,                 # process this top hits
	'id|idfile:s'   => \$idfile,
	'uh|use_hash' 		 => \$use_hash,
	'nt|notimer' => \$notimer,
	'nc|nocheck' => \$no_check,
	'lines:i'    => \$hash_lines,
	'x|maskxn:s'		=> \$maskxn,
	'extract:s'	=> \$extract,	# file to extract fasta from into CDS
	'single' => \$single_copy
);
if ( !@ARGV ) {
	warn("\nUnless you give me at least one BLAST file, there is nothing for me to do!\n");
	pod2usage;

}
if ($use_hash) { undef($store_hash); }
my @blastfiles = @ARGV;
if ($debug2) { $debug = 1; }

if ($extract){
	unless (-s $extract){die ("Could not find $extract");}
	unless (-f "$extract.index"){
		print "Indexing $extract...\n";
        	my $inx = Bio::Index::Fasta->new(-filename => "$extract.index",-write_flag => 1);
	        $inx->make_index($extract);
        }
        $extract_inx = Bio::Index::Fasta->new(-filename => "$extract.index") || die ("Could not get index for $extract\n");
}


if ($idfile && -s $idfile) {
        my $pattern;
        $pattern='/^\s*(\S+)\s+/';
        my @test_lines=`head $idfile`;
        foreach my $test (@test_lines){if ($test=~/^>/){$pattern="Bio::SeqIO";}}
        print "Building hash from $idfile with $pattern\n";

        if ($pattern eq "Bio::SeqIO"){
                my $id_obj=new Bio::SeqIO( -file => $idfile,-format => "fasta" );
                while ( my $seq = $id_obj->next_seq() ) {
                        $ids{$seq->id()}=1;
                }
        }
        else{
                open( IN, $idfile ) || die();
                my $flag = 0;
                while ( my $line = <IN> ) {
                        if ( $line =~ $pattern ) {
                                $ids{$1} = 1;
                                if ( $flag == 0 ) {
                                        print "Hash presence of $idfile verified\n";
                                        $flag = 1;
                                }
                        }
                }
              close(IN);
        }
}

foreach my $blastfile (@blastfiles) {
	if ($use_hash && $blastfile=~/^(.+)\.hash$/){$blastfile=$1;}
	unless ( -s $blastfile || -s "$blastfile.hash" ) {
		die("Can't find $blastfile\n");
	}
	print "Processing $blastfile...\n";
	my $logfile=$blastfile.".analysis.".$cut_score.".".$cut_evalue;
	#reset global varialbes
	$query_tlength           = int(0);
	$query_tlength_with_hits = int(0);
	if ($use_hash) {
		unless ( -e "$blastfile.hash" ) {
			die("You requested to use a previous hash but I can't find $blastfile.hash\n" );}
	}
	if ( $debug || $debug2 ) {
		use Data::Dumper;
		$debugfile = $blastfile . ".debug";
		open( DEBUG, ">$debugfile" );
	}
	open( LOG, ">$logfile" );
	if ($store_hash) {
		if ( -s "$blastfile.hash" ) {
			warn "$blastfile.hash already exists. Are you sure you want to rebuild the HASH-table?\nWait 1sec if yes; Ctl-C otherwise\n";
			sleep(1);
		}
		open( HASH, ">$blastfile.hash" );
		print HASH "TYPE\tSEQ NAME\tSEQ LENGTH\tGLOBAL HSP No\tIDENTICAL\tCONSERVED\tSTART\tEND\tBITSCORE\tEVALUE\tLENGTH PROP\tLOCAL HSP No\tLOCAL HIT No\tDIRECTION\n";
	}
	&process_blast($blastfile);
	if ($store_hash) { close(HASH); }

	if ($extract){
		print "Preparing CDS and UTR files...\n";
		my $fasta_cds = new Bio::SeqIO( -file => ">$extract.putative_cds.$blastfile", -format => "fasta" );
		$fasta_cds->width(15000);
		my $fasta_utr = new Bio::SeqIO( -file => ">$extract.putative_utr.$blastfile", -format => "fasta" );
		$fasta_utr->width(15000);
		foreach my $query_name (keys %extr_hash){
			my $gene_obj = $extract_inx->fetch($query_name);
			die "Cannot find sequence $query_name in $extract\n" if !$gene_obj;
			die "Sequence seems to be a protein (or have IUPAC codes), cannot create CDS/UTR files!\n" if ($gene_obj->alphabet() ne 'dna') ;
			if (!$gene_obj){next;}
			my $gene_length=$gene_obj->length();
			my $id=$gene_obj->id();
			my $start=$extr_hash{$query_name}{"start"};
			my $end=$extr_hash{$query_name}{"end"};
			my $direction = $extr_hash{$query_name}{"direction"};
			#$gene_obj = $gene_obj->revcom() if $direction eq 'R';
			if (!$start){$start=$gene_length+1;}
			if (!$end){$end=0;}
			#print "Processing id $id length $gene_length. CDS start $start end $end \n";
			my $cds=new Bio::Seq();
			my $utr5=new Bio::Seq();
			my $utr3=new Bio::Seq();
			my $cds_seq = $direction eq 'R' ? &revcomp($gene_obj->subseq($start,$end)) : $gene_obj->subseq($start,$end);
			$cds->seq($cds_seq);
			my $new_id=$id;
			my $desc = "CDS ".$start."_".$end;
			$desc.=' R' if $direction eq 'R';
			$cds->desc($desc);
			$cds->id($new_id);
			$fasta_cds->write_seq($cds);
			if ($start>2){
				my $new_start=1;
				my $new_end=$start-1;
				my $seq = $direction eq 'R' ? &revcomp($gene_obj->subseq($new_start,$new_end)) : $gene_obj->subseq($new_start,$new_end);
				$utr5->seq($gene_obj->subseq($new_start,$new_end));
				my $new_id=$id;
				my $desc = "5UTR ".$new_start."_".$new_end.' F' if $direction eq 'F';
				$desc = "3UTR ".$new_start."_".$new_end.' R' if $direction eq 'R';
				$utr5->desc($desc);
				$utr5->id($new_id);
				$fasta_utr->write_seq($utr5);
				$fasta_utr->write_seq($utr5);
			}
			if ($end<$gene_length){
				my $new_start=$end+1;
				my $new_end=$gene_length;
				$utr3->seq($gene_obj->subseq($new_start,$new_end));
				my $new_id=$id;
				my $desc = "3UTR ".$new_start."_".$new_end.' F' if $direction eq 'F';
				$desc = "5UTR ".$new_start."_".$new_end.' R' if $direction eq 'R';
				$utr5->desc($desc);
				$utr3->id($new_id);
				$fasta_utr->write_seq($utr3);
			}
		}
	}



	#print $print_statement;
	#print LOG $print_statement;
	close(LOG);
	close(DEBUG);
	print "All Done! See $logfile\n#==============================#\n\n";

}
##################################################################
# SUBROUTINES
##################################################################
sub process_blast($) {
	my $blastfile = shift;
	print "Parsing BLAST report $blastfile\n";
	unless ($use_hash) {
		my $total_queries = `grep -c Query= $blastfile`;
		chomp($total_queries);
		#$timer->attr( min => 0, max => $total_queries );
		#$timer->restart;

		# verify blast has run to completion
		unless ($no_check) {
			if ( $report_style eq 'blast' || $report_style eq 'BLAST' ) {
				my $result;
				my @lines = `tail -n 20 $blastfile`;
				foreach my $ln (@lines){
					$result = 1 if $ln=~/^Matrix/;
				}
				if ( !$result ) {
					print("Sorry, but it seems your BLAST report is incomplete\n");
					return;
				}
			}
		}
	}
	my ( $hash_ref_db_elements, $hash_ref_queries, $db_name, $db_length,
		 $db_entries, $query_total );
	my $annotation_redundancy=int(0);
# to reduce memory we can split the read hash in two rounds and build each hash independently then emptying it.
# but we only want to read the blast file once so the read blast does no longer build the hash, but it does print it
# to be read later
	if ($use_hash) {
		print "Reading database hash\n";
		(
		   $hash_ref_db_elements, $db_name, $db_length, $db_entries,
		   $query_total
		) = &read_hash( $blastfile, "database" );
	} else {
		print "Building new database hash\n";
		(
		   $hash_ref_db_elements, $db_name, $db_length, $db_entries,
		   $query_total
		) = &build_hash($blastfile);
	}
	
	
	## Now with hashes built continue to process things
	my ( $db_ided, $query_ided );

	# prepare arrays for the Stats
	my $db_with_hit=int(0);
	my ($conserved_array_ref, $identical_array_ref, $unique_matches_ref,$all_matches_ref,$length_proportions_ref);
	if ($hash_ref_db_elements) {
		print "Processing database data\n";
		my $print_statement;
		# memory explodes here.
		(	 $conserved_array_ref, $identical_array_ref, $unique_matches_ref,
			 $all_matches_ref,  $length_proportions_ref,   $db_with_hit
		) = parse_blasthash( $hash_ref_db_elements, "database",$blastfile );
		$db_ided = $db_with_hit;
		my $db_with_hit_ratio =	  sprintf( "%d %%", $db_with_hit / $db_entries * 100 );
		undef(%$hash_ref_db_elements);    # is this correct way to empty memory?
		undef($hash_ref_db_elements);
		$print_statement .="\nYou searched versus\t$db_name:\nTotal entries:\t$db_entries\nDatabase positions:\t$db_length positions\nIdentified:\t$db_with_hit\t($db_with_hit_ratio).\n";
		print ".";
		foreach my $prop (sort keys %{$length_proportions_ref}){
			$print_statement .="Database entries identified with total length at least \t".$prop."%\t".$length_proportions_ref->{$prop}."\n";
		}
		print ".";
		my @unique_matches = @$unique_matches_ref;
		my ( $unique_sum, $unique_mean, $unique_sd, $unique_median ) =  &prepare_stats( \@unique_matches, "unique" );
		$unique_mean = sprintf( "%d", $unique_mean );
		my $prop_db_identified =
		  sprintf( "%d %%", $unique_sum / $db_length * 100 );
		$print_statement .="Non-overlapping positions identified:\t$unique_sum ($prop_db_identified), with mean length $unique_mean (SD=$unique_sd) and median $unique_median\n";
		undef(@unique_matches);
		print ".";
		my @all_matches = @$all_matches_ref;
		my ( $all_sum, $all_mean, $all_sd, $all_median ) = &prepare_stats( \@all_matches, "all" );
		$all_mean = sprintf( "%d", $all_mean );
		my $unique_ratio = sprintf( "%.6f", $all_sum / $unique_sum );
		$annotation_redundancy+=$unique_ratio;
		$print_statement .="Overlapping positions identified:\t$all_sum, with mean length $all_mean (SD=$all_sd) and median $all_median.\nOverlapping/non-overlapping ratio:\t$unique_ratio\n";
		undef(@all_matches);
		print ".";
		my @conserved_array = @$conserved_array_ref;
		my ( $conserved_sum, $conserved_mean, $conserved_sd, $conserved_median ) = &prepare_stats( \@conserved_array, "conserved" );
		$conserved_mean = sprintf( "%.2f %%", $conserved_mean );
		$print_statement .="Mean of conserved positions:\t$conserved_mean (SD=$conserved_sd), and median $conserved_median\n";
		undef(@conserved_array);
		print ".";
		my @identical_array = @$identical_array_ref;
		my ( $identical_sum, $identical_mean, $identical_sd, $identical_median )
		  = &prepare_stats( \@identical_array, "identical" );
		$identical_mean = sprintf( "%.2f %%", $identical_mean );
		$print_statement .="Mean of identical positions:\t$identical_mean (SD=$identical_sd), and median $identical_median\n";
		undef(@identical_array);
		print ". Done\n";
		if ($debug) { &create_debug_log( $hash_ref_db_elements, "database" ); }
		print LOG $print_statement;
	}
	( $hash_ref_queries, $db_name, $db_length, $db_entries, $query_total ) =  &read_hash( $blastfile, "queries" );
	  
	#QUERY data
	# update Query length if X masking requested
	if ($maskxn){
		$query_tlength=int(0);
		if (int($maskxn)){
			print "Query length provided by user as $maskxn.\n";
			$query_tlength=$maskxn;
		}
		else {
			my $fasta_obj=new Bio::SeqIO( -file => $maskxn,-format => "fasta" ) || die("$maskxn is not a fasta file\n");
			while ( my $seq = $fasta_obj->next_seq() ) {
				my $sequence=$seq->seq();
				$sequence=~s/X+//ig;
				$query_tlength+=length($sequence);
			}
			print "Query length estimated from inputfile as query_tlength. No X/Ns\n";
		}
		if (!$query_tlength || $query_tlength==0){
			die "Query length is 0. This should not have happened. Did you give the -maskx argument properly?\n";
		}
	}
# \@conserved_array, \@identical_array, \@unique_matches,\@all_matches,  \%length_proportions,   $element_number
	if ($hash_ref_queries) {
		print "Processing query data\n";
		my $print_statement;
		my ( $unique_matches_ref, $all_matches_ref, $length_proportions_ref ,$queries_with_hit) = parse_blasthash( $hash_ref_queries, "queries",$blastfile );
		$query_ided = $queries_with_hit;
		my $queries_with_hit_ratio =  sprintf( "%d %%", $queries_with_hit / $query_total * 100 );
		undef(%$hash_ref_queries);
		undef($hash_ref_queries);
		$print_statement .="\nFrom the query perspective:\nTotal queries: $query_total. Queries with hit:\t$queries_with_hit ($queries_with_hit_ratio)\n";
		print '.';
        foreach my $prop (sort keys %{$length_proportions_ref}){
            $print_statement .="Query sequences that are covered by reference at least \t".$prop."%\t".$length_proportions_ref->{$prop}."\n";
        }
        print '.';
		#my @conserved_array=@$conserved_array_ref;
		#my @identical_array=@$identical_array_ref;
		my @unique_matches = @$unique_matches_ref;
		my @all_matches    = @$all_matches_ref;
		my ( $unique_sum, $unique_mean, $unique_sd, $unique_median ) =	  &prepare_stats( \@unique_matches, "unique" );
		$unique_mean = sprintf( "%d", $unique_mean );
		my $prop_query_identified = sprintf( "%d %%", $unique_sum / $query_tlength * 100 );
		$print_statement .= "Total length:\t$query_tlength positions.\nNon-overlapping positions identified:\t$unique_sum ($prop_query_identified), with mean length $unique_mean (SD=$unique_sd) and median $unique_median\n";
		undef(@unique_matches);
		print ".";
		my ( $all_sum, $all_mean, $all_sd, $all_median ) =	  &prepare_stats( \@all_matches, "all" );
		$all_mean = sprintf( "%d", $all_mean );
		my $unique_ratio = sprintf( "%.6f ", $all_sum / $unique_sum );
		$annotation_redundancy+=$unique_ratio;
		$print_statement .="Overlapping positions identified:\t$all_sum, with mean length $all_mean (SD=$all_sd) and median $all_median.\nOverlapping/non-overlapping ratio:\t$unique_ratio\n";
		undef(@all_matches);
		print ". Done\n";
		$print_statement .="\nTotal annotation redundancy:\t$annotation_redundancy\n";
		if ($debug) { &create_debug_log( $hash_ref_queries, "Queries" ); }
		print LOG $print_statement;
	}
	return (0);    # success!
}

sub prepare_stats($$) {
	my $array_ref  = shift;
	my $array_type = shift;
	my @array      = @$array_ref;
	my $print_statement;
	if ($debug) { print DEBUG "$array_type\n"; print DEBUG Dumper @array; }
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@array);
	my $mean = $stat->mean();
	if ($mean) { $mean = sprintf( "%.4f", $mean ); }
	else       { $mean = 0; }
	my $sd = $stat->standard_deviation();
	if ($sd) { $sd = sprintf( "%.2f", $sd ); }
	else     { $sd = 0; }
	my $median = $stat->median();
	if ( !$median ) { $median = 0; }
	my $sum = $stat->sum();
	return ( $sum, $mean, $sd, $median );
}

sub parse_blasthash ($$$) {
	my $hash_ref  = shift;   # this either the db elements or the query elements
	my $hash_type = shift;   # database or queries
	my $blastfile = shift;
	my @conserved_array;     # for database only
	my @identical_array;     # for database only
	my @unique_matches;
	my @all_matches;
	my %length_proportions;
	my $element_number = int(0);    # hit or queries present in hash
	open (FULL_LENGTH,">$blastfile.$hash_type.full");
	foreach my $element ( keys %$hash_ref ) {

		# this is not really necessary, unless we want to build a graph later on (or debug)
		if ($debug2) {
			my $length = $hash_ref->{$element}{"length"};
			for ( my $i = 1 ; $i <= $length ; $i++ ) {
				$hash_ref->{$element}{"pos"}{$i} = int(0);
			}
		}
		# for each database element, we have top aln_prop
		#unless ( $hash_type eq "queries" ) {
			if ($hash_ref->{$element}{"aln_prop"}){
				print FULL_LENGTH $element."\t".$hash_ref->{$element}{'align_name'}."\t".$hash_ref->{$element}{"aln_prop"}."\n" if $hash_ref->{$element}{"aln_prop"}>=0.80;
				$length_proportions{25}++ if $hash_ref->{$element}{"aln_prop"}>=0.25;
				$length_proportions{50}++ if $hash_ref->{$element}{"aln_prop"}>=0.50;
				$length_proportions{60}++ if $hash_ref->{$element}{"aln_prop"}>=0.60;
				$length_proportions{70}++ if $hash_ref->{$element}{"aln_prop"}>=0.70;
				$length_proportions{75}++ if $hash_ref->{$element}{"aln_prop"}>=0.75;
				$length_proportions{80}++ if $hash_ref->{$element}{"aln_prop"}>=0.80;
				$length_proportions{90}++ if $hash_ref->{$element}{"aln_prop"}>=0.90;
				$length_proportions{95}++ if $hash_ref->{$element}{"aln_prop"}>=0.95;
			}
		#}
		# memory explode
		# go to each hsp and get data out

		foreach my $hsp ( keys %{ $hash_ref->{$element}{"hsp"} } ) {
			# first do the positions of element.
			my $element_start_position =  $hash_ref->{$element}{"hsp"}{$hsp}{"start"};
			my $element_end_position =  $hash_ref->{$element}{"hsp"}{$hsp}{"end"};
			for ( my $i = $element_start_position ;
				  $i <= $element_end_position ;
				  $i++ ){
				my $old_size = $hash_ref->{$element}{"pos"}{$i} if ($debug2);
				$hash_ref->{$element}{"pos"}{$i}++;
				if ($debug2) {
					my $new_size = $hash_ref->{$element}{"pos"}{$i};
					print DEBUG "\nFound one! $element ($element_start_position,$element_end_position): $i was $old_size and is $new_size\n";
				}
			}

			# now build up conserved/identical arrays
			unless ( $hash_type eq "queries" ) {
				push( @conserved_array,
					  $hash_ref->{$element}{"hsp"}{$hsp}{"conserved"} );
				push( @identical_array,
					  $hash_ref->{$element}{"hsp"}{$hsp}{"identical"} );
			}
		}

# we finished cycling through the HSPs of this element so if we want to do stats to the positions of this element
# only, then this is the place to do it. At the moment we don't; so move on.
	}
	close (FULL_LENGTH);
	print "Hash parsed. Populating arrays...\n";

# Positions on subjects have been stored so we can now accumulate the results and populate the arrays
# Now with positions, we could see how the exact coverage per base is too... for future implementation;
	foreach my $element ( keys %$hash_ref ) {
		$element_number++;
		foreach my $position ( keys %{ $hash_ref->{$element}{"pos"} } ) {
			$hash_ref->{$element}{"aggregate_total"} +=
			  $hash_ref->{$element}{"pos"}{$position};
			$hash_ref->{$element}{"unique_total"}++;
		}

		# avoid undefs going in!
		my $match  = $hash_ref->{$element}{"aggregate_total"};
		my $unique = $hash_ref->{$element}{"unique_total"};
		if ($match)  { push( @all_matches,    $match ); }
		if ($unique) { push( @unique_matches, $unique ); }
	}
	if ( $hash_type eq "queries" ) {
		return ( \@unique_matches, \@all_matches, \%length_proportions,$element_number );
	} elsif ( $hash_type eq "database" ) {
		return (
				 \@conserved_array, \@identical_array, \@unique_matches,
				 \@all_matches,  \%length_proportions,   $element_number
		);
	}
}

# use stored hash
sub read_hash($$) {
	my $blastfile = shift;
	my $hash_type = shift;
	warn "Reading hash $hash_type\n";
	my ( %hash, $db_name, $db_length, $db_entries, $query_total );
	print "Stored HASH found, parsing...\n";
	open( HASH, "$blastfile.hash" );
	my (%hit_counting, %query_counting);
	my $total_hit_counter=int(0);
	my $hsp_counter=int(0);
	my $query_counter=int(0);
	#bookmark
	# timer
	my $line_num=int(0);
	if ($hash_lines) { $line_num = $hash_lines; }
	elsif (-s "$blastfile.hash") {
		$line_num = `wc $blastfile.hash`;
		chomp($line_num);
		$line_num =~ s/^\s*(\d+)\s.+$/$1/;
	}
	print "HASH has $line_num lines...\n";
	my $line_counter;
	#$timer->attr( min => 0, max => $line_num );
	#$timer->restart;
	LINE: while ( my $line = <HASH> ) {
		$line_counter++;
		#if ( !$notimer && $line_counter =~ /0000$/ ) {
		#	print $timer->report( "eta: %E min, %40b %p\r", $line_counter );
		#}
		if ( $line =~ /^\#/ ) {
			if ( $line =~ /^\#Total.+\:(\d+)/ ) {
				$query_tlength_with_hits = $1;
			} elsif ( $line =~ /^\#Overall.+\:(\d+)/ ) {
				$query_tlength = $1;
			} elsif ( $line =~ /^\#Queries total:(\d+)/ ) {
				$query_total = $1;
			} elsif ( $line =~ /^\#DB/ ) {
				if    ( $line =~ /name:(\S+)/ )     { $db_name    = $1; }
				elsif ( $line =~ /length\:(\d+)/ )  { $db_length  = $1; }
				elsif ( $line =~ /entries\:(\d+)/ ) { $db_entries = $1; }
			} else {
				next;
			}
		}
		chomp($line);
		my @data = split( "\t", $line );
		if ( @data && $hash_type eq "database" && $data[0] eq "HIT" ) {
			my $db_element_name   = $data[1];
			my $db_element_length = $data[2];
			$hsp_counter                 = $data[3];
			my $frac_id_db        = $data[4];
			my $frac_cons_db      = $data[5];
			my $hstart            = $data[6];
			my $hend              = $data[7];
			my $score             = $data[8];
			my $eval              = $data[9];
			my $ref_aln_prop      = $data[10];
			my $local_hsp_counter = $data[11];
			my $local_hit_counter = $data[12] if $data[12];
			my $direction = $data[13];
			my $query_name = $data[14];
			next if ($limit && $local_hit_counter && $limit < $local_hit_counter);
			# a) to prevent low scoring HSPs to be accepted - repeats; b) to allow reparsing of hash. Next line.
			if ( $score < $cut_score || $eval > $cut_evalue ) {next;} 
			
			$hsp_counter++;
			$hit_counting{$db_element_name} =
			  1;    #just so we can get a value of how many hits we have...
			$hash{$db_element_name}{"length"} = $db_element_length;
			$hash{$db_element_name}{"hsp"}{$hsp_counter} = {
												   "identical" => $frac_id_db,
												   "conserved" => $frac_cons_db,
												   "start"     => $hstart,
												   "end"       => $hend,
												   "score"     => $score,
												   "evalue"    => $eval,
			};
			$hash{$db_element_name}{"aln_prop"} = $ref_aln_prop if ($local_hsp_counter==1 &&(!$hash{$db_element_name}{"aln_prop"} || $hash{$db_element_name}{"aln_prop"}<$ref_aln_prop) );
			$hash{$db_element_name}{'align_name'} =  $query_name;
		} elsif ( @data && $hash_type eq "queries" && $data[0] eq "QUERY" ) {
			my $query_name   = $data[1];
			my $query_length = $data[2];
			$hsp_counter  = $data[3];
			my $qstart = $data[6];
			my $qend   = $data[7];
			my $score  = $data[8];
			my $eval   = $data[9];
			my $query_aln_prop       = $data[10];
			my $local_hsp_counter = $data[11];
			my $local_hit_counter = $data[12] if $data[12];
			my $direction = $data[13];
			my $hit_name = $data[14];
			next if ($limit && $local_hit_counter && $limit < $local_hit_counter);
			if ( $score < $cut_score ) {
				next;
			} # a) to prevent low scoring HSPs to be accepted - repeats; b) to allow reparsing of hash
			if ( $eval > $cut_evalue ) { next; }
			$query_counting{$query_name} =
			  1;    #just so we can get a value of how many hits we have...
			$hash{$query_name}{"length"} = $query_length;

			if ($extract){
				my $start = $qstart;
				my $end = $qend;
				if ($end<$start){my $t=$end;$end=$start;$start=$t;}
				$extr_hash{$query_name}{"start"}=$start if !$extr_hash{$query_name}{"start"} || $start < $extr_hash{$query_name}{"start"};
				$extr_hash{$query_name}{"end"}=$end if !$extr_hash{$query_name}{"end"} || $end > $extr_hash{$query_name}{"end"};
				$extr_hash{$query_name}{"direction"}= $data[13];
#die "query $query_name has end ".$extr_hash{$query_name}{"end"};
			}

			$hash{$query_name}{"hsp"}{$hsp_counter} = {
				"start"  => $qstart,
				"end"    => $qend,
				"score"  => $score,
				"evalue" => $eval,
			};
			$hash{$query_name}{"aln_prop"} = $query_aln_prop if ($local_hsp_counter==1 &&(!$hash{$query_name}{"aln_prop"} || $hash{$query_name}{"aln_prop"}<$query_aln_prop) );
			$hash{$query_name}{'align_name'} =  $hit_name;
		}
	}
	close(HASH);
	foreach my $key ( keys %hit_counting )   { $total_hit_counter++; }
	foreach my $key ( keys %query_counting ) { $query_counter++; }
	#my $elapsed = $timer->report("%L");
	#print "\nTime elapsed: $elapsed min.\n";
	#print LOG "\nTime elapsed: $elapsed min.\n";
	#print "Using cut-off evalue of $cut_evalue and bit-score $cut_score.\nCalculating Stats...\n";
	#print LOG "Using cut-off evalue of $cut_evalue and bit-score $cut_score.\n";
	my $hash_ref = \%hash;
	return ( $hash_ref, $db_name, $db_length, $db_entries, $query_total );
}

# build new hash
sub build_hash ($) {
	my $blastfile = shift;
	my (
		 %hash_db_elements,        %hash_queries,      $db_name,
		 $query_counter,           $total_hit_counter, 
		 $query_tlength_with_hits, $db_length,         $db_entries,
		 $query_total,             $global_hit_index
	);
	my $hitless = int(0);
	my $blast_obj =	  Bio::SearchIO->new( -file => $blastfile, -format => $report_style );
	my $hsp_counter=int(0); #global HSP counter
	print "Building HASH...\n";
	while ( my $result = $blast_obj->next_result() ) {

		# this is the timer
		$query_total++;
		#if ( !$notimer && $query_total =~ /0000$/ ) {
		#	print $timer->report( "eta: %E min, %40b %p\r", $query_total );
		#}
		if ($idfile){
			unless ( exists $ids{$result->query_name} ){
				next;
			}
		}
		my $query_length = $result->query_length();
		$query_tlength += $query_length;
		my $hitcount = $result->num_hits;
		if ( $hitcount == 0 ) { $hitless++; next; } # skip queries with no hits.

		# shall we allow multiple blast reports in one file? Don't think so....
		unless ($db_name) {
			$db_name = $result->database_name();
			$db_name =~ s/\s+$//;
		}
		unless ($db_length) {
			$db_length = $result->database_letters();
			$db_length =~ s/\D//g;
		}
		unless ($db_entries) {
			$db_entries = $result->database_entries();
			$db_entries =~ s/\D//g;
		}
		$query_tlength_with_hits += $query_length;
		my $query_name = $result->query_name();
		$hash_queries{$query_name}{"length"} = $query_length;
		my $hit_counter       = int(0);
		my $query_hsp_counter = int(0);

# db_element is, essentially, each element in the BLAST database. Contrast with Query which is the elements in the query dataset.		
		while ( my $hit = $result->next_hit ) {
			$hit_counter++;    # for this query
			last if ( $limit && $limit < $hit_counter );
			$hit->overlap(5);
			my ($qcontigs, $scontigs) = Bio::Search::SearchUtils::tile_hsps($hit);
			my $reference_length = $hit->length();
			my $tscore = $hit->bits();
			my $teval  = $hit->significance();
			if ( $tscore < $cut_score ) { next; }
			if ( $teval > $cut_evalue ) { next; }
			if ( $hit->rank() == 1 ) { $query_counter++; }
			$global_hit_index++;
			my $db_element_name   = $hit->name();
			my $db_element_length = $hit->length();
			my $qstrand = $hit->strand('query');
			my $rstrand = $hit->strand('hit');
			my $direction = $qstrand == $rstrand ? 'F' : 'R';
			my ($reference_aln_length,$query_aln_length);
			$hash_db_elements{$db_element_name}{"length"} = $db_element_length;
			if ($scontigs==1){
				my $hsp = $hit->hsp()||next;
				my $hsp_rank = $hsp ->rank();
				$reference_aln_length=$hsp->length('hit');
				my $ref_aln_prop = sprintf("%.4f",$reference_aln_length/$reference_length);
				$query_aln_length=$hsp->length('query');
				my $query_aln_prop = sprintf("%.4f",$query_aln_length/$query_length);
				my $score = $hsp->bits();
				my $eval  = $hsp->evalue();
				#if ( $hsp_rank == 1 && $score < $cut_score ) { next; }
				#if ( $hsp_rank == 1 && $eval > $cut_evalue ) { next; }
				$hsp_counter++;        # global hsps index in whole blast report
				$query_hsp_counter++;  # Index for query HSPs
				my ( $hstart, $hend ) = $hsp->range('hit');
				my $frac_id_db   = $hsp->frac_identical('hsp');
				my $frac_cons_db = $hsp->frac_conserved('hsp');
				$frac_id_db   = sprintf( "%.4f", $frac_id_db );
				$frac_cons_db = sprintf( "%.4f", $frac_cons_db );
				$hash_db_elements{$db_element_name}{"hsp"}{$hsp_counter} = {
					"identical" => $frac_id_db,
					"conserved" => $frac_cons_db,
					"start"     => $hstart,
					"end"       => $hend,
					"score"     => $score,
					"evalue"    => $eval
				};
				$hash_queries{$query_name}{'align_name'} = $db_element_name;
				$hash_queries{$query_name}{"aln_prop"} = $query_aln_prop if (
					$hsp_rank==1 && (!$hash_queries{$query_name}{"aln_prop"} || $hash_queries{$query_name}{"aln_prop"}<$query_aln_prop)
				 );
				$hash_db_elements{$db_element_name}{"aln_prop"} = $ref_aln_prop if ($hsp_rank==1 &&(!$hash_db_elements{$db_element_name}{"aln_prop"} || $hash_db_elements{$db_element_name}{"aln_prop"}<$ref_aln_prop) );
				$hash_db_elements{$db_element_name}{'align_name'} = $query_name;
				my ( $qstart, $qend ) = $hsp->range('query');
				my $frac_id_query =   "N/A";    # actually we dont want these values for queries
				my $frac_cons_query = "N/A";
				if ($store_hash) {
					print HASH "HIT\t$db_element_name\t$db_element_length\t$hsp_counter\t$frac_id_db\t$frac_cons_db\t$hstart\t$hend\t$score\t$eval\t$ref_aln_prop\t$hsp_rank\t$hit_counter\tF\t$query_name\n";
					print HASH "QUERY\t$query_name\t$query_length\t$hsp_counter\t$frac_id_query\t$frac_cons_query\t$qstart\t$qend\t$score\t$eval\t$query_aln_prop\t$hsp_rank\t$hit_counter\t$direction\t$db_element_name\n";
				}
			}else{
				foreach my $contig (@{$scontigs}){
					$reference_aln_length+= abs($contig->{'stop'}-$contig->{'start'})+1;
					my $ref_aln_prop = sprintf("%.4f",$reference_aln_length/$reference_length);
					my $query_aln_prop = sprintf("%.4f",$reference_aln_length/$query_length);
					my $hsp=@{$contig->{'hsps'}}[0];
					my $hsp_rank = $hsp ->rank();
					my $score = $hsp->bits();
					my $eval  = $hsp->evalue();
					#if ($hsp_rank ==1 && $score < $cut_score ) { next; }
					#if ($hsp_rank ==1 && $eval > $cut_evalue ) { next; }
					$hsp_counter++;        # global hsps index in whole blast report
					$query_hsp_counter++;  # Index for query HSPs
					my ( $hstart, $hend ) = $hsp->range('hit');
					my $frac_id_db   = $hsp->frac_identical('hsp');
					my $frac_cons_db = $hsp->frac_conserved('hsp');
					$frac_id_db   = sprintf( "%.4f", $frac_id_db );
					$frac_cons_db = sprintf( "%.4f", $frac_cons_db );
					$hash_db_elements{$db_element_name}{"hsp"}{$hsp_counter} = {
						"identical" => $frac_id_db,
						"conserved" => $frac_cons_db,
						"start"     => $hstart,
						"end"       => $hend,
						"score"     => $score,
						"evalue"    => $eval
					};
					$hash_queries{$query_name}{'align_name'} = $db_element_name;
					$hash_queries{$query_name}{"aln_prop"} = $query_aln_prop if ($hsp_rank==1 &&(!$hash_queries{$query_name}{"aln_prop"} || $hash_queries{$query_name}{"aln_prop"}<$query_aln_prop) );
					$hash_db_elements{$db_element_name}{"aln_prop"} = $ref_aln_prop if ($hsp_rank==1 &&(!$hash_db_elements{$db_element_name}{"aln_prop"} || $hash_db_elements{$db_element_name}{"aln_prop"}<$ref_aln_prop) );
					$hash_db_elements{$db_element_name}{'align_name'} = $query_name;
					my ( $qstart, $qend ) = $hsp->range('query');
					my $frac_id_query =   "N/A";    # actually we dont want these values for queries
					my $frac_cons_query = "N/A";
					if ($store_hash) {
						print HASH "HIT\t$db_element_name\t$db_element_length\t$hsp_counter\t$frac_id_db\t$frac_cons_db\t$hstart\t$hend\t$score\t$eval\t$ref_aln_prop\t$hsp_rank\t$hit_counter\tF\t$query_name\n";
						print HASH "QUERY\t$query_name\t$query_length\t$hsp_counter\t$frac_id_query\t$frac_cons_query\t$qstart\t$qend\t$score\t$eval\t$query_aln_prop\t$hsp_rank\t$hit_counter\t$direction\t$db_element_name\n";
					}
				}
			}
		}
	}
	if ($store_hash) {
		print HASH
		  "#Total length of queries with hits:$query_tlength_with_hits\n";
		print HASH "#Overall query length:$query_tlength\n";
		print HASH "#Queries total:$query_total\n";
		print HASH
		  "#DB name:$db_name\n#DB length:$db_length\n#DB entries:$db_entries\n";
	}
	foreach my $key ( keys %hash_db_elements ) {
		$total_hit_counter++;
	}
	#my $elapsed = $timer->report("%L");
	#print
#"\nTime elapsed: $elapsed min.\nFound $query_counter queries with $total_hit_counter unique hits. $global_hit_index non-unique hits totalling $hsp_counter HSPs using cut-off evalue of $cut_evalue and bit-score $cut_score.\n$hitless queries had no hits and have been discarded.\nCalculating Stats...\n";
	#print LOG
#"\nTime elapsed: $elapsed min.\nFound $query_counter queries with $total_hit_counter unique hits. $global_hit_index non-unique hits totalling $hsp_counter HSPs using cut-off evalue of $cut_evalue and bit-score $cut_score.\n$hitless queries had no hits and have been discarded.\n";
	my $hash_ref_db_elements = \%hash_db_elements;
	my $hash_ref_queries     = \%hash_queries;

	# no longer return the query hash, empty query info after printing it out.
	return ( $hash_ref_db_elements, $db_name, $db_length, $db_entries,
			 $query_total );
}

sub create_debug_log ($$) {
	my $hash_ref  = shift;
	my $hash_type = shift;
	print DEBUG "\n\nLooking at $hash_type\n";
	foreach my $element ( keys %$hash_ref ) {
		print DEBUG "\nelement $element\t";
		foreach my $position ( sort { $a <=> $b }
							   ( keys %{ $hash_ref->{$element}{"pos"} } ) )
		{
			my $counter = $hash_ref->{$element}{"pos"}{$position};
			print DEBUG "\n\t$position has $counter";
		}
	}
}


sub revcomp {
  my $dna = shift;
  my $revcomp = reverse(uc($dna));
  $revcomp =~ tr/ACGT/TGCA/;
  return $revcomp;
}

