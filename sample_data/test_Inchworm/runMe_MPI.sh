#!/bin/sh
if [ ! -e jellyfish.kmers.fa ]; then
    gunzip -c jellyfish.kmers.fa.gz > jellyfish.kmers.fa
fi

mpirun -n 2 ../../Inchworm/src/MPIinchworm  --kmers jellyfish.kmers.fa --token tok  -K 25 --DS > mpi.iworm.out
