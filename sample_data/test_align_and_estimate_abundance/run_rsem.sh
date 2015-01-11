if [ ! -e Trinity.fasta ]; then
    gunzip -c Trinity.fasta.gz > Trinity.fasta
fi

if [ ! -e reads.left.fq ]; then
    gunzip -c reads.left.fq.gz > reads.left.fq
fi

if [ ! -e reads.right.fq ]; then
    gunzip -c reads.right.fq.gz > reads.right.fq
fi


../../util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads.left.fq --right reads.right.fq --gene_trans_map Trinity.fasta.gene_trans_map --aln_method bowtie --est_method RSEM --output_dir rsem_estimate --SS_lib_type RF --prep_reference
