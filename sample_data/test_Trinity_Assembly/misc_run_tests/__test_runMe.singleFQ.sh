#!/bin/bash -ve

## Error: ./__test_runMe.singleFQ.sh: line 44: ../../util/RSEM_util/summarize_RSEM_fpkm.pl: No such file or directory

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

../../Trinity.pl --seqType fq --single reads.left.fq.gz  --JM 1G  --CPU 4  --output trinity_single_outdir



sleep 2

######################################################
## align reads back to the transcripts using Bowtie ##
######################################################

sleep 2

../../util/alignReads.pl --single reads.left.fq --target trinity_single_outdir/Trinity.fasta --aligner bowtie --seqType fq -o bowtie_single_outdir

##### Done aligning reads #######

sleep 2

###########################################
# use RSEM to estimate read abundance  ####
###########################################

sleep 2

../../util/RSEM_util/run_RSEM.pl --transcripts trinity_single_outdir/Trinity.fasta --name_sorted_bam bowtie_single_outdir/bowtie_single_outdir.nameSorted.bam

###### Done running RSEM ########

sleep 2


############################
# compute FPKM values   ####
############################

sleep 2

../../util/RSEM_util/summarize_RSEM_fpkm.pl --transcripts trinity_single_outdir/Trinity.fasta --RSEM RSEM.isoforms.results --fragment_length 76 --group_by_component | tee Trinity.RSEM.fpkm

#############################
####   Done.  ###############
#############################


