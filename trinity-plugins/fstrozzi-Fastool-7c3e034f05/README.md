Fastool
=======

A simple and quick tool to read huge FastQ and FastA files (both normal and gzipped) and manipulate them.

It makes use of the KSeq library (http://lh3lh3.users.sourceforge.net/kseq.shtml) for fast access to FastQ/A files.

Installation
------------

Clone this repository and run make.

Usage
-----

     Usage: %s (--rev) (--append [string_to_append_to_header]) (--to-fasta) (--illumina-trinity) sequences_1.fastq/a sequences_2.fastq/a ...

--to-fasta (optional): convert FastQ files to FastA format.

--rev (optional): reverse complement all the sequences in the dataset (both FastQ and FastA).

--append (optional): add a string at the end of each sequence header (both FastQ and FastA).

--illumina-trinity (optional): directly converts Casava 1.8+ FastQ ID format to Trinity Fasta input format (appending /1 and /2 for PE reads)

Examples
--------

FastQ conversion for Trinity pipeline
    
    fastool --illumina-trinity sequences.fastq > sequences.fasta

FastQ to FastA conversion

    fastool --to-fasta sequences.fastq > sequences.fasta

Return the reverse complement

    fastool --rev sequences.fastq > reverse_complement.fastq

Append '/1' to the end of the sequence ID

    fastool --append /1 sequences.fasta > forward_sequences.fasta

Append '/2' to the end of the sequence ID and return the reverse complement

    fastool --append /2 --rev sequences.fasta > reverse_sequences.fasta

Can read compressed files directly

    fastool --to-fasta sequences.fastq.gz > sequences.fasta

Can process more then one file

    fastool --to-fasta --append /1 sequences1.fastq sequences2.fastq sequences3.fastq > all_sequences.fasta
    
Can be used with pipes

    cat sequences.fastq | fastool --to-fasta --rev > reverse_complement.fasta        

License
-------

     Permission is hereby granted, free of charge, to any person obtaining
     a copy of this software and associated documentation files (the
     "Software"), to deal in the Software without restriction, including
     without limitation the rights to use, copy, modify, merge, publish,
     distribute, sublicense, and/or sell copies of the Software, and to
     permit persons to whom the Software is furnished to do so, subject to
     the following conditions:

     The above copyright notice and this permission notice shall be
     included in all copies or substantial portions of the Software.

     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
     EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
     MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
     NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
     BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
     ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
     CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
     SOFTWARE.

Copyright
---------

Copyright (c) 2012 Francesco Strozzi

