#!/bin/bash
# Prepare reads (first 10M pairs) and reference transcript sets for the FL regression benchmark.
# usage: prep_data.sh <data_dir>
# Requires: gsutil (reads are in a requester-readable public bucket), curl, python3, gzip, bzip2.
set -eu
D=${1:?usage: prep_data.sh <data_dir>}
HERE=$(cd $(dirname $0) && pwd)
mkdir -p $D/reads $D/ref
cd $D/reads

B=gs://ctat-public-ftp-resources/trinity_test_rnaseq_datasets
N=40000000   # 10M fastq records = 40M lines

# head closes the pipe early, so the upstream exit status is 141 (SIGPIPE); that's expected.
[ -s mouse_10M.left.fq  ] || { gsutil cat $B/mouse.all.Left.fq.gz  | gzip  -dc | head -n $N > mouse_10M.left.fq  || true; }
[ -s mouse_10M.right.fq ] || { gsutil cat $B/mouse.all.Right.fq.gz | gzip  -dc | head -n $N > mouse_10M.right.fq || true; }
[ -s pombe_10M.left.fq  ] || { gsutil cat $B/SP2.all.Left.fq.bz2   | bzip2 -dc | head -n $N > pombe_10M.left.fq  || true; }
[ -s pombe_10M.right.fq ] || { gsutil cat $B/SP2.all.Right.fq.bz2  | bzip2 -dc | head -n $N > pombe_10M.right.fq || true; }
for f in *_10M.*.fq; do n=$(wc -l < $f); [ $n -eq $N ] || { echo "ERROR: $f has $n lines, expected $N" >&2; exit 1; }; done

cd $D/ref
# Pinned references (record versions in results_history.tsv if these ever change)
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M39
ENSFUNGI=https://ftp.ensemblgenomes.ebi.ac.uk/pub/fungi/release-63/fasta/schizosaccharomyces_pombe/cdna
POMBASE=https://www.pombase.org/data/genome_sequence_and_features/feature_sequences

[ -s gencode.vM39.transcripts.fa.gz ]    || curl -sSfL -o gencode.vM39.transcripts.fa.gz    $GENCODE/gencode.vM39.transcripts.fa.gz
[ -s gencode.vM39.pc_transcripts.fa.gz ] || curl -sSfL -o gencode.vM39.pc_transcripts.fa.gz $GENCODE/gencode.vM39.pc_transcripts.fa.gz
[ -s pombe.cdna.all.fa.gz ]              || curl -sSfL -o pombe.cdna.all.fa.gz              $ENSFUNGI/Schizosaccharomyces_pombe.ASM294v2.cdna.all.fa.gz
[ -s pombe.cds.fa.gz ]                   || curl -sSfL -o pombe.cds.fa.gz                   $POMBASE/cds.fa.gz   # not versioned by PomBase; check md5 in README

# FL_trans_analysis_pipeline.pl requires headers of the form  >transcript_id;gene_id
# mouse cDNA: protein-coding transcripts (GENCODE header field 8 = transcript_type)
zcat gencode.vM39.transcripts.fa.gz | awk -F'|' '/^>/{keep=($8=="protein_coding"); if(keep) print $1";"$2; next} keep{print}' > mouse.cdna.ref.fa
# mouse CDS: CDS span of each pc transcript, identical CDSs collapsed per gene
python3 -I $HERE/make_cds_ref.py gencode.vM39.pc_transcripts.fa.gz mouse.cds.ref.fa
# pombe cDNA (Ensembl Fungi, includes annotated UTRs)
zcat pombe.cdna.all.fa.gz | awk '/^>/{t=substr($1,2); g=$4; sub(/^gene:/,"",g); print ">" t ";" g; next} {print}' > pombe.cdna.ref.fa
# pombe CDS (PomBase), protein-coding only (drop SPNCRNA/snRNA/snoRNA/rRNA/tRNA entries)
zcat pombe.cds.fa.gz | awk '/^>/{t=substr($1,2); keep=(t !~ /^SP(NCRNA|SNRNA|SNORNA|RRNA|TRNA)/); if(keep) print ">" t ";" t; next} keep{print}' > pombe.cds.ref.fa

md5sum *.gz
grep -c '>' *.ref.fa
