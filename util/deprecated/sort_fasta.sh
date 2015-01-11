#!/bin/bash


if [ ! $1 ]
then
    echo
    echo usage $0 file.fasta \> file.sorted.fasta
    echo
    exit 1
fi


cat $1 | tr '\n>' '$\n'  | sort | tr '\n$' '>\n' | sed  '$d' 


exit $?
