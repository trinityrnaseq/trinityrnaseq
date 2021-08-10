
VERSION=`cat VERSION.txt`

docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq:$VERSION Trinity  --seqType fq --left `pwd`/test_data/reads.left.fq.gz  --right `pwd`/test_data/reads.right.fq.gz  --max_memory 1G --CPU 4 --output `pwd`/trinity_out_dir
