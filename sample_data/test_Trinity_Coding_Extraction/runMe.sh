#!/bin/sh 

if [ -e Trinity.fasta.gz ] && [ ! -e Trinity.fasta ]
then
    gunzip -c Trinity.fasta.gz > Trinity.fasta
fi


../../trinity-plugins/transdecoder/TransDecoder -t Trinity.fasta -m 50

echo
echo 
echo See best_candidates.\*  for candidate ORFs
echo



