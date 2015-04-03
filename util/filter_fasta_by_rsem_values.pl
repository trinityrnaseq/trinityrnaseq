#!/usr/bin/env perl

=head1 NAME

filter_fasta_by_rsem_values.pl - Use RSEM relative abundance values to filter a
transcript assembly FASTA file

=head1 SYNOPSIS

USAGE: filter_fasta_by_rsem_values.pl 
            --rsem_output=/path/to/RSEM.isoforms.results[,...]
            --fasta=/path/to/Trinity.fasta
            --output=/path/to/output.fasta
          [ --tpm_cutoff=1.0
            --fpkm_cutoff=0.5
            --isopct_cutoff=0.05
          ]

=head1 OPTIONS

B<--rsem_output,-r>
    This file is the RSEM output file, usually from util/RSEM_util/run_RSEM.pl
    Provide a comma-separated list for multiple RSEM output files, or specify
    parameter like so:
    -r RSEM.isoforms.A -r RSEM.isoforms.B -r RSEM.isoforms.C ...

B<--fasta,-f>
    The FASTA file representing transcript assemblies from the primary Trinity run.

B<--output,-o>
    The output FASTA file to be created. These are the sequences whose scores are
    greater or equal to the user-specified cutoffs.

B<--filtered_output,-e>
    Optional. Pass this if you want a FASTA file created of those sequences which
    do NOT meet the user-specified cutoffs.

B<--tpm_cutoff,-t>
    Optional.  Will filter transcripts, keeping those with TPM values equal to or
    greater than this.  See the INPUT section for more detail here.

B<--fpkm_cutoff,-c>
    Optional.  Will filter transcripts, keeping those with FPKM values equal to or
    greater than this.  See the INPUT section for more detail here.

B<--isopct_cutoff,-i>
    Optional.  Will filter transcripts, keeping those with isoform percentage values
    equal to or greater than this.  See the INPUT section for more detail here.

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

Trinity contains a pipeline for read alignment, visualization and relative abundance 
estimation here:

    http://trinityrnaseq.sourceforge.net/analysis/align_visualize_quantify.html

Among the products of this is the file 'RSEM.isoforms.results', which contains calculated
values such as TPM, FPKM and IsoPct.  See INPUT for more on this.

The user can use this RSEM file to filter the source transcript assembly FASTA file
based on any combination of these cutoffs.

=head1  INPUT

The input RSEM file looks like this:

    transcript_id   gene_id length  effective_length        expected_count  TPM     FPKM    IsoPct
    comp100579_c0_seq1      comp100579_c0   204     10.71   0.00    0.00    0.00    0.00
    comp100585_c0_seq1      comp100585_c0   301     85.75   0.00    0.00    0.00    0.00
    comp100592_c0_seq1      comp100592_c0   439     222.65  9.00    0.39    0.31    100.00
    comp100599_c0_seq1      comp100599_c0   268     54.93   0.00    0.00    0.00    0.00
    comp100606_c0_seq1      comp100606_c0   807     590.56  36.65   0.60    0.47    16.24
    comp100606_c0_seq2      comp100606_c0   783     566.56  181.35  3.10    2.44    83.76

The corresponding input FASTA has headers like:

    >comp100579_c0_seq1 len=204 path=[1:0-203]
    >comp100585_c0_seq1 len=301 path=[1:0-88 90:89-300]
    >comp100592_c0_seq1 len=439 path=[1:0-438]
    >comp100599_c0_seq1 len=268 path=[1:0-267]
    >comp100606_c0_seq1 len=807 path=[1:0-285 784:286-309 287:310-806]
    >comp100606_c0_seq2 len=783 path=[1:0-285 287:286-782]

To do the actual filtering, you need to pass at least one of the cutoff parameters.
They are additive - if more than one are passed, each will be applied to futher filter
the input set.

There are filters for TPM, FPKM, and PctIso.  From the RSEM documentation, they are:

'TPM' stands for Transcripts Per Million. It is a relative measure of transcript abundance. The
sum of all transcripts' TPM is 1 million. 'FPKM' stands for Fragments Per Kilobase of transcript
per Million mapped reads. It is another relative measure of transcript abundance.  'IsoPct' stands
for isoform percentage. It is the percentage of this transcript's abandunce over its parent gene's
abandunce. If its parent gene has only one isoform or the gene information is not provided, this
field will be set to 100.


=head1  OUTPUT

The output is a FASTA file with the same headers and sequence, only filtered based
on the parameters passed.

=head1  CONTACT

    Joshua Orvis
    jorvis@gmail.com

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
                          'rsem_output|r=s@',
                          'fasta|f=s',
                          'output|o=s',
                          'tpm_cutoff|t=s',
                          'fpkm_cutoff|c=s',
                          'isopct_cutoff|i=s',
                          'filtered_output|e=s',
                          'log|l=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

if (@ARGV) {
    die "Error, did not recognize parameters: @ARGV ";
}

## make sure everything passed was peachy
&check_parameters(\%options);

## open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">$options{log}") || die "can't create log file: $!";
}

_log( "INFO: cutoffs: fpkm_cutoff=(" . ($options{fpkm_cutoff} || '') . "), perc_comp_fpkm_cutoff=(" . ($options{perc_comp_fpkm_cutoff} || '') . ")");

_log( "INFO: Opening RSEM file ($options{rsem_output})" );
my $rsem = load_rsem_output( $options{rsem_output} );

_log( "INFO: creating output file: ($options{output})" );
open(my $ofh, ">$options{output}") || logdie("ERROR: failed to create output file: $!");
open (my $rsem_ofh, ">$options{output}.rsem") or logdie $!;
print $rsem_ofh join("\t", "transcript_id", "gene_id", "length", "effective_length",
                     "expected_count",  "TPM", "FPKM", "IsoPct", "iso_per_gene", "rsem_filename") . "\n"; # header line

my $filtered_ofh;
if ( $options{filtered_output} ) {
    _log( "INFO: creating optional filtered output file: ($options{filtered_output})" );
    open($filtered_ofh, ">$options{filtered_output}") || logdie("ERROR: failed to create filtered output file: $!");
}

_log( "INFO: reading input FASTA file: ($options{fasta})" );
open( my $ifh, $options{fasta} ) || logdie("ERROR: failed to read input FASTA file: $!");

my $keep = 0;


while ( my $line = <$ifh> ) {
    if ( $line =~ /^\>(\S+)/ ) {
        my $seq_id = $1;
        
        ## make sure we have cutoffs for this
        if ( ! exists $$rsem{$seq_id} ) {
            logdie("ERROR: found a transcript ID ($seq_id) in the FASTA that wasn't in the RSEM file");
        }
        
        my $gene_id = $$rsem{iso_to_gene}->{$seq_id} or logdie("ERROR, cannot determine gene_id for $seq_id");
        my $num_iso_per_gene = $$rsem{iso_count_per_gene}->{$gene_id} or logdie("ERROR, cannot determine number of isoforms for gene: $gene_id");
        
        $keep = 1; ## optimism first.  just keep swimming
        
        
        my @rsem_entries = @{$$rsem{$seq_id}}; # this is why people hate or love perl.
        
        my @rsem_entries_meet_requirements;

        foreach my $rsem_entry (@rsem_entries) {
            ## if any individual rsem entry for this transcript meets all specified filtering criteria, then we keep it.
            my $local_keep = 1;

            $local_keep = 0 if ( defined $options{fpkm_cutoff} && $rsem_entry->{fpkm} < $options{fpkm_cutoff} );
            $local_keep = 0 if ( defined $options{tpm_cutoff} && $rsem_entry->{tpm} < $options{tpm_cutoff} );
            $local_keep = 0 if ( defined $options{isopct_cutoff} && $rsem_entry->{isopct} < $options{isopct_cutoff} 
                                 && $num_iso_per_gene > 1);
            
            if ($local_keep) {
                push (@rsem_entries_meet_requirements, $rsem_entry);
            }
        }

        if (@rsem_entries_meet_requirements) {
            
            foreach my $rsem_entry (@rsem_entries_meet_requirements) {
                print $rsem_ofh $rsem_entry->{line} . "\t$num_iso_per_gene\t" . $rsem_entry->{file} . "\n";
            }
            
        }
        else {
            $keep = 0;
        }
        
    }
    
    print $ofh $line if $keep;
    if ( $keep == 0 && $options{filtered_output} ) {
        print $filtered_ofh $line;
    }
}

close $ofh;
close $rsem_ofh;

exit(0);


sub load_rsem_output {
    my $files_aref = shift;
    
    my @files = split(/,/,join(',', @$files_aref));
    
    ## relative abundance data.  looks like:
    #   $rel{'comp3119_c0_seq1'} = { fpkm   => 63.55,
    #                                tpm    => 10324,
    #                                isopct => 31.43
    #                              }
    my %rel = ();
    
    foreach my $file (@files) {
        
        open(my $ifh, $file) || logdie("ERROR: failed to read rsem_output file: $!");
        
        while (my $line = <$ifh>) {
            chomp $line;
            
            next if $line =~ /^\s*$/;
            my @cols = split("\t", $line);
            
            if ( scalar @cols == 8 && $cols[3] ne 'effective length' ) {
                if ( exists $rel{$cols[0]} && grep { $_->{file} eq $file } @{$rel{$cols[0]}} ) { 
                    logdie("ERROR: found more than one entry in the RSEM file for $cols[0]");
                } else {
                    my $struct = { fpkm => $cols[6], 
                                   tpm => $cols[5], 
                                   isopct => $cols[7],
                                   line => $line,
                                   file => $file,
                                   trans => $cols[0],
                                   gene => $cols[1], 
                               };
                    
                    push (@{$rel{$cols[0]}}, $struct);
                }
                
                
            }
        }
            
    }

    _log("INFO: Loaded RSEM values for " . scalar(keys %rel) . " transcripts");


    my %trans_to_gene;
    my %gene_to_iso;
    
    ## generate trans to gene mapping and count isoforms per gene.
    foreach my $trans_id (keys %rel) {
        my @rsem_entries = @{$rel{$trans_id}};
        

        foreach my $rsem_entry (@rsem_entries) {
            my $gene_id = $rsem_entry->{gene};
            $gene_to_iso{$gene_id}->{$trans_id} = 1;
            $trans_to_gene{$trans_id} = $gene_id;
        }
    }

    $rel{iso_to_gene} = \%trans_to_gene;
    
    foreach my $gene_id (keys %gene_to_iso) {
        my @trans = keys %{$gene_to_iso{$gene_id}};
        my $num_trans = scalar(@trans);
        $rel{iso_count_per_gene}->{$gene_id} = $num_trans;
    }

    return \%rel;
}


sub _log {
    my $msg = shift;

    print $logfh "$msg\n" if $logfh;
}

sub logdie {
    my $msg = shift;
    
    print STDERR "$msg\n";
    print $logfh "$msg\n" if $logfh;
    
    exit(1);
}


sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( rsem_output fasta output );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            logdie("ERROR: --$option is a required option");
        }
    }
    
    ## at least one of the cutoffs must be defined
    if ( ! defined $$options{fpkm_cutoff} &&
         ! defined $$options{tpm_cutoff} &&
         ! defined $$options{isopct_cutoff}) {
        logdie("ERROR: You must define at least one cutoff for filtering.  See the PERLDOC");
    }
}
