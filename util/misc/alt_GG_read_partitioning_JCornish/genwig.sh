#!/bin/bash
sam2wig(){
    chrlen=$1
    chr=$2
    infi=$3
    nosingle=$4
    minins=$5
    maxins=$6
    threads=$7

    outfi="$infi.$chr.wig"
    tmp="$infi.$chr.sam"

    
    samtools view -@ $threads $infi $chr | \
    ./genwig2.py $chrlen $chr - $outfi $nosingle $minins $maxins
    echo "completed: $chr"
}
export -f sam2wig

while getopts ":b:m:n:o:p:q:" opt; do
    case $opt in
    b)
        if [ ! -f $OPTARG ]
            then
                echo "ERROR: Unable to find input bam $OPTARG" >&2
                exit 1
            else
                inbam="$OPTARG"
        fi
    ;;
    m)
        maxins=$OPTARG
    ;;
    n)
        minins=$OPTARG
    ;;
    o)
        if [ ! -d $OPTARG ];
            then
                echo "ERROR: Unable to find output directory $OPTARG" >&2
                exit 1
            else
                outdir="$OPTARG"
        fi
    ;;
    p)
        nproc=$OPTARG
    ;;
    q)
        nosingle=$OPTARG
    ;;
    \?)
        echo "ERROR: Unknown argument $opt" >&2
        exit 1
    ;;
    esac
done

#dump chromosome info from bam index
chrinfo="$inbam.chr"
#sort is required to preserve order between trinity and here
samtools idxstats $inbam | grep -v "^*" | cut -f1,2 | sort > $chrinfo

#generate wigs for segs/chrs
parallel -j $nproc --xapply \
sam2wig {1} {2} $inbam $nosingle $minins $maxins $nproc \
::: `cut -f2 $chrinfo` ::: `cut -f1 $chrinfo`

#cat to final wig and remove file
if [ -f "$inbam.wig" ]
then
    rm "$inbam.wig"
fi

while read chr len
do
    cat "$inbam.$chr.wig" >> "$inbam.wig"
    rm "$inbam.$chr.wig"
done < $chrinfo
