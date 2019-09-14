singularity exec -e /seq/RNASEQ/TOOLS/TRINITY/SINGULARITY/Trinity.simg  Trinity --seqType fq --single `pwd`/reads.left.fq.gz  --max_memory 1G --CPU 2 --output `pwd`/trinity_out_dir_singularity

