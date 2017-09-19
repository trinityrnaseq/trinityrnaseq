#!/usr/bin/env python

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


def run_cmd(cmd):
    logger.info("Running: " + " ".join(cmd))
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (output, error) = process.communicate()
    logger.debug(output)
    if process.returncode != 0:
        logger.error(error)
        exit("Error while running command \"" + str(cmd) + "\":\n" + error)
        exit("Error while running command \"" + " ".join(cmd) + "\":\n" + error)


def main():

    FORMAT = "%(asctime)-15s %(levelname)s %(module)s.%(name)s.%(funcName)s at %(lineno)d :\n\t%(message)s\n"
    global logger 
    logger = logging.getLogger()
    logging.basicConfig(filename='variant_calling.log', format=FORMAT, filemode='w', level=logging.DEBUG)
    # add a new Handler to print all INFO and above messages to stdout
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    logger.addHandler(ch)
    logger.info('Start')

    parser = argparse.ArgumentParser(description="This script requires you have the following dependencies:\nSamtools: \"samtools\" in your path\nJava: \"java\" in your path\nPicard-Tools: \"$PICARD_HOME\" with the path to Picard-Tools's bin\nSTAR: \"STAR\" in your path\nGATK: \"$GATK\" with the path to GATK's bin\n", epilog="", formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-st', '--supertranscript', dest="supertranscript_file", type=str, required=True, help="Path to the SuperTranscripts fasta file.")

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-p', '--paired', dest="paired_reads", type=str, nargs=2, help="Pair of paired ends read files.")

    group.add_argument('-s', '--single', dest="single_reads", type=str, help="Single reads file.")

    parser.add_argument('-o', '--output', dest="out_path", type=str, required=True, help="Path to the folder where to generate the output.")

    parser.add_argument('-l', '--readlength', dest="readlength", type=int, help="Size of the reads (used for STAR --sjdbOverhang). If no value is specified the length of the first read will be used.")

    parser.add_argument('-t', '--threads', dest="nthreads", type=str, default="4", help="Number of threads to use for tools that are multithreaded.")

    parser.add_argument('-m', '--maxram', dest="maxram", type=str, default="50000000000", help="Maximum amount of RAM allowed for STAR's genome generation step (only change if you get an error from STAR complaining about this value).")

    args = parser.parse_args()

    # if args.single_reads:
    #     if args.left_reads or args.right_reads:
    #         exit("Specify either a --single read file or a pair of --left and --right read files.")
    # elif not args.left_reads or not args.right_reads:
    #     exit("Please specify a correct input for read files.")

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

    st_path = os.path.realpath(args.supertranscript_file)

    if not args.readlength:
        with open(reads_paths[0], "r") as f:
            f.readline()  # skip header
            sjdbOverhang = str(len(f.readline().rstrip()) - 1)
    else:
        sjdbOverhang = str(args.readlength - 1)

    # check if output directory exists, if not create
    real_path = os.path.realpath(args.out_path)
    if not os.path.isdir(real_path):
        os.makedirs(real_path)

    # generate supertranscript index
    os.chdir(os.path.dirname(st_path))
    logger.info("Generating SuperTranscript index.")
    cmd = shlex.split("samtools faidx " + os.path.basename(st_path))
    run_cmd(cmd)

    # generate supertranscript Picard dictionary
    logger.info("Generating Picard dictionary.")
    if not ".".join(os.path.basename(st_path).split(".")[:-1]) + ".dict" in os.listdir("."):
        cmd = shlex.split("java -jar " + PICARD_HOME + "/picard.jar CreateSequenceDictionary R= " + os.path.basename(st_path) + " O= " + ".".join(os.path.basename(st_path).split(".")[:-1]) + ".dict")
        run_cmd(cmd)

    # move to output folder
    os.chdir(real_path)

    # generate genome folder for STAR's first pass
    if not "star_genome_generate" in os.listdir("."):
        logger.info("Generating genome folder for STAR")
        if not os.path.isdir("star_genome"):
            os.mkdir("star_genome")
        os.chdir("star_genome")
        cmd = shlex.split("STAR --runThreadN " + args.nthreads + " --runMode genomeGenerate  --genomeDir ./  --genomeFastaFiles " + st_path + " --limitGenomeGenerateRAM " + args.maxram)
        run_cmd(cmd)
        open("../star_genome_generate", 'w').close()
        os.chdir("../")
    else:
        logger.info("STAR's first pass genome folder already generated, skipping.")

    # run STAR's alignment
    if not "star_alignment_run" in os.listdir("."):
        logger.info("Running STAR's first pass.")
        if not os.path.isdir("star_alignment"):
            os.mkdir("star_alignment")
        os.chdir("star_alignment")
        cmd = str("STAR --runThreadN " + args.nthreads
                  + " --genomeDir ../star_genome  "
                  + " --twopassMode Basic "
                  + " --readFilesIn " + " ".join(reads_paths))
        
        if re.search("\.gz$", reads_paths[0]):
            cmd += " --readFilesCommand zcat "

        cmd = shlex.split(cmd)
        run_cmd(cmd)
        open("../star_alignment_run", 'w').close()
        os.chdir("../")
    else:
        logger.info("STAR's first pass already processed, skipping.")


    # clean and convert sam file with Picard-Tools
    logger.info("Cleaning and Converting sam file with Picard-Tools.")
    if not os.path.isdir("picard_steps"):
        os.mkdir("picard_steps")
    os.chdir("picard_steps")

    if not "picard_addorreplace" in os.listdir(".."):
        cmd = shlex.split("java -jar " + PICARD_HOME + " AddOrReplaceReadGroups I=../star_pass_2/Aligned.out.sam O=rg_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample")
        run_cmd(cmd)
        open("../picard_addorreplace", 'a').close()
    else:
        logger.info("Picard's AddOrReplaceReadGroups already processed, skipping.")        

    if not "picard_markduplicates" in os.listdir(".."):
        cmd = shlex.split("java -jar " + PICARD_HOME + "/picard.jar MarkDuplicates I=rg_added_sorted.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics")
        run_cmd(cmd)
        open("../picard_markduplicates", 'a').close()
    else:
        logger.info("Picard's MarkDuplicates already processed, skipping.")    

    if not os.path.isdir("../gatk_steps"):
        os.mkdir("../gatk_steps")
    os.chdir("../gatk_steps")

    if not "gatk_splitncigar" in os.listdir(".."):
        cmd = shlex.split("java -jar " + os.getenv("GATK_HOME") + "/GenomeAnalysisTK.jar -T SplitNCigarReads -R " + st_path + " -I ../picard_steps/dedupped.bam -o split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS")
        run_cmd(cmd)
        open("../gatk_splitncigar", 'a').close()
    else:
        logger.info("Picard's SplitNCigarReads already processed, skipping.")    

    # do the actual variant calling
    logger.info("Variant Calling using Haplotype Caller.")
    if not "gatk_haplotypecaller" in os.listdir(".."):
        cmd = shlex.split("java -jar " + GATK_HOME + "/GenomeAnalysisTK.jar -T HaplotypeCaller -R " + st_path + " -I ./split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o output.vcf")
        run_cmd(cmd)
        open("../gatk_haplotypecaller", 'a').close()
    else:
        logger.info("GATK's HaplotypeCaller already processed, skipping.")    

    # do some basic filtering
    logger.info("Doing some basic filtering of vcf.")
    if not "gatk_variantfiltration" in os.listdir(".."):
        cmd = shlex.split("java -jar " + GATK_HOME + "/GenomeAnalysisTK.jar -T VariantFiltration -R " + st_path + " -V output.vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o filtered_output.vcf")
        run_cmd(cmd)
        open("../gatk_variantfiltration", 'a').close()
    else:
        logger.info("GATK's VariantFiltration already processed, skipping.")    
    logger.info("Done!")


if __name__ == "__main__":
    main()
