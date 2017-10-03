#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import os
import sys
import argparse
import subprocess
import shlex
import logging
import re

sys.path.insert(0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)),
                                     "..", "..", "..", "PyLib"]))

import Pipeliner


logger = None


def main():

    FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
    global logger 
    logger = logging.getLogger()
    logging.basicConfig(filename='variant_calling.log', format=FORMAT, filemode='w', level=logging.DEBUG)
    # add a new Handler to print all INFO and above messages to stdout
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)

    parser = argparse.ArgumentParser(
        description=str("This script requires you have the following dependencies:\n" +
                        "Samtools: \"samtools\" in your path\n" +
                        "Java: \"java\" in your path\n" +
                        "Picard-Tools: env var \"$PICARD_HOME\" with the path to Picard-Tools's bin\n" +
                        "STAR: \"STAR\" in your path\n" +
                        "GATK: env var \"$GATK_HOME\" with the path to GATK's bin\n"),
        epilog="", formatter_class=argparse.RawTextHelpFormatter) 
    
    parser.add_argument('--st_fa', '--supertranscript_fasta', dest="st_fa", type=str, required=True, help="Path to the SuperTranscripts fasta file.")

    parser.add_argument('--st_gtf', '--supertranscript_gtf', dest="st_gtf", type=str, required=True, help="Path to the SuperTranscript gtf file.")

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-p', '--paired', dest="paired_reads", type=str, nargs=2, help="Pair of paired ends read files.")

    group.add_argument('-s', '--single', dest="single_reads", type=str, help="Single reads file.")

    parser.add_argument('-o', '--output', dest="out_path", type=str, required=True, help="Path to the folder where to generate the output.")

    parser.add_argument('-l', '--sjdbOverhang', dest="sjdbOverhang", default=150, type=int, help="Size of the reads (used for STAR --sjdbOverhang). default=150")
    
    parser.add_argument('-t', '--threads', dest="nthreads", type=str, default="4", help="Number of threads to use for tools that are multithreaded.")

    parser.add_argument('-m', '--maxram', dest="maxram", type=str, default="50000000000", help="Maximum amount of RAM allowed for STAR's genome generation step (only change if you get an error from STAR complaining about this value).")

    args = parser.parse_args()

    

    PICARD_HOME = os.getenv("PICARD_HOME")
    if not PICARD_HOME:
        exit("Error, missing path to Picard-Tools in $PICARD_HOME.")

    GATK_HOME = os.getenv("GATK_HOME")
    if not GATK_HOME:
        exit("Error, missing path to GATK in $GATK.")



    # get real paths before changing working directory in case they are relative paths
    if args.paired_reads:
        reads_paths = [os.path.realpath(f) for f in args.paired_reads]
    else:
        reads_paths = [os.path.realpath(args.single_reads)]

    st_fa_path = os.path.realpath(args.st_fa)

    st_gtf_path = os.path.realpath(args.st_gtf)
    
    
    # check if output directory exists, if not create
    real_path = os.path.realpath(args.out_path)
    if not os.path.isdir(real_path):
        os.makedirs(real_path)

    # move to output folder
    os.chdir(real_path)
    
    checkpoint_dir = os.path.abspath(os.path.basename(st_fa_path)) + ".gatk_chkpts"
    pipeliner = Pipeliner.Pipeliner(checkpoint_dir)
            
    # generate supertranscript index
    logger.info("Generating SuperTranscript index.")
    pipeliner.add_commands([Pipeliner.Command("samtools faidx {}".format(st_fa_path), "samtools_faidx_st.ok")])
    pipeliner.run()
    
    # generate supertranscript Picard dictionary
    logger.info("Generating Picard dictionary.")
    dict_file = re.sub("\.[^\.]+$", ".dict", st_fa_path)
    if os.path.isfile(dict_file):
        open(checkpoint_dir + "/picard_dict_st.ok", 'a').close()
    else:
        pipeliner.add_commands([Pipeliner.Command("java -jar " + PICARD_HOME + "/picard.jar" +
                                                  " CreateSequenceDictionary R=" + st_fa_path +
                                                  " O=" + dict_file,
                                                  "picard_dict_st.ok")])
        pipeliner.run()

    # generate genome folder for STAR's first pass
    logger.info("Generating genome folder for STAR")
    
    star_genome_generate_cmd = str("STAR --runThreadN " +
                                   args.nthreads +
                                   " --runMode genomeGenerate" +
                                   " --genomeDir star_genome_idx " +
                                   " --genomeFastaFiles {} ".format(st_fa_path) +
                                   " --sjdbGTFfile {} ".format(st_gtf_path) +
                                   " --sjdbOverhang {} ".format(args.sjdbOverhang) +
                                   " --limitGenomeGenerateRAM {}".format(args.maxram) )
        
    pipeliner.add_commands([
        Pipeliner.Command("mkdir star_genome_idx", "mkdir_star_genome_idx.ok"),

        Pipeliner.Command(star_genome_generate_cmd,
                          "star_genome_generate.ok")
        ])
    pipeliner.run()

    
    # run STAR's alignment
    logger.info("Running STAR alignment.")
    cmd = str("STAR --runThreadN " + args.nthreads
              + " --genomeDir star_genome_idx "
              + " --runMode alignReads "
              + " --twopassMode Basic "
              + " --alignSJDBoverhangMin 10 "
              + " --outSAMtype BAM SortedByCoordinate "
              + " --limitBAMsortRAM {} ".format(args.maxram) 
              + " --readFilesIn " + " ".join(reads_paths) )
    
    if re.search("\.gz$", reads_paths[0]):
        cmd += " --readFilesCommand 'gunzip -c' "
        
    pipeliner.add_commands([Pipeliner.Command(cmd, "star_aln.ok")])
    pipeliner.run()

    
    
    # clean and convert sam file with Picard-Tools
    logger.info("Cleaning and Converting sam file with Picard-Tools.")

    pipeliner.add_commands([
        Pipeliner.Command("java -jar " + PICARD_HOME + "/picard.jar " +
                          " AddOrReplaceReadGroups " +
                          "I=Aligned.sortedByCoord.out.bam " +
                          "O=rg_added_sorted.bam " +
                          "SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample",
                          "add_read_groups.ok"),

        Pipeliner.Command("java -jar " + PICARD_HOME + "/picard.jar " +
                          " MarkDuplicates " +
                          "I=rg_added_sorted.bam O=dedupped.bam " +
                          "CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics",
                          "mark_dups.ok"),

        Pipeliner.Command("java -jar " + GATK_HOME + "/GenomeAnalysisTK.jar " +
                          "-T SplitNCigarReads -R " + st_fa_path +
                          " -I dedupped.bam -o splitNCigar.bam " +
                          " -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS",
                          "splitNCigarReads.ok")
        ])
    pipeliner.run()

    
    
    # do the actual variant calling
    logger.info("Variant Calling using Haplotype Caller.")

    pipeliner.add_commands([
        Pipeliner.Command("java -jar " + GATK_HOME + "/GenomeAnalysisTK.jar " +
                          "-T HaplotypeCaller -R " + st_fa_path +
                          " -I ./splitNCigar.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf",
                          "haplotypecaller.ok")
        ])
    pipeliner.run()
    
    # do some basic filtering
    logger.info("Doing some basic filtering of vcf.")
    
    pipeliner.add_commands([
        Pipeliner.Command("java -jar " + GATK_HOME + "/GenomeAnalysisTK.jar " +
                          "-T VariantFiltration -R " + st_fa_path +
                          " -V output.vcf -window 35 -cluster 3 " +
                          "-filterName FS -filter \"FS > 30.0\" " +
                          "-filterName QD -filter \"QD < 2.0\" -o filtered_output.vcf",
                          "variant_filt.ok")
        ])

    pipeliner.run()

    
    logger.info("Done!")


if __name__ == "__main__":
    main()
