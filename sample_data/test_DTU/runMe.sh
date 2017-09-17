../../Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl --genes_fasta data/minigenome.fa --genes_gtf data/minigenome.gtf --samples_file data/samples.txt --out_prefix G

rm -f ./*.ok

../../Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py  --trinity_fasta data/mini.Trinity_fmt.fasta --out_prefix trinSuper

../../Analysis/SuperTranscripts/DTU/dexseq_wrapper.pl --genes_fasta trinSuper.fasta --genes_gtf trinSuper.gtf --samples_file data/samples.txt --out_prefix S

./compare_dexseq_results.pl G.results.dat S.results.dat > compare.dat
