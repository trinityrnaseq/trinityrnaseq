#!/usr/bin/env python

#    use_express.py
#    October 2012  Matthew MacManes (macmanes@gmail.com)
#
#   This wrapper is free software: you can redistribute it and/or modify
#
#  v.0.3.1 Changes: Do not remake index if already made, added -k30 option to Bowtie2 mapping step


import sys
import subprocess
import optparse
import shutil
import os
from datetime import datetime, date, time
from Bio import SeqIO
import os.path

print ""
print ""
print ""
print "******************************************************************"
print "***   run_express.py v0.3.1                                 ******"
print "***   To run this program, you must have bowtie2 and eXpress******"
print "***   installed and in your $PATH                           ******"
print "******************************************************************"
print ""



##########################################
## date function
##########################################

def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")

##########################################
## Options
##########################################
def getOptions():
    parser = optparse.OptionParser(usage="usage: python %prog -b input.fa  -t index_name -p [num threads] -X [insert size] -l left.fq -r right.fq -n file.sam]",
                          version="%prog 0.3.1")
    parser.add_option("-b", "--b2base",
                      dest="b2base",
                      default="Trinity.fasta",
                      metavar='file.fa',
                      help="fasta file for B2 index (?Trinity.fasta)")
    parser.add_option("-t", "--target",
                      dest="target",
                      metavar='index',
                      default="index",
                      help="Name of bowtie2 index",)
    parser.add_option("-p", "--threads",
                      dest="threads",
                      metavar='[INT]',
                      default="2",
                      help="Number of threads to use",)
    parser.add_option("-X", "--insert",
                      dest="insert",
                      default="500",
                     metavar='[INT]', 
                     help="Max inner distance",)
    parser.add_option("-l", "--left",
                      dest="left",
                      metavar='file.fq',
                      default="",
                      help="comma sep list of left reads",)
    parser.add_option("-r", "--right",
                      dest="right",
                      metavar='file.fq',
                      default="",
                      help="comma sep list of right reads",)
    parser.add_option("-o", "--outdir",
                      dest="outdir",
                      metavar='path to output directory',
		      default=".",
                      help="output directory",)
    parser.add_option("-n", "--name",
                      dest="name",
                      metavar='SAM filename',
		      default="hits.sam",
                      help="full path filename for SAM file from bowtie2")
    parser.add_option("-U", "--unpaired",
                      dest="unpaired",
                      metavar='unpaired reads',
		      default="",
                      help="full path to unpaired reads")
    (options, args) = parser.parse_args()

    return options


##########################################
## alignment procedure
##########################################
#def numbering(options, awker):
#	with open('%s.tmp' %(options.b2base),'w') as stdout_fh:
#		num = subprocess.Popen(['awk', awker, options.b2base], stdout=stdout_fh)
#		output = num.communicate()
def b2build(options):
	b2b = subprocess.Popen(['bowtie2-build', '%s' % (options.b2base), options.target], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	output = b2b.communicate()
	assert b2b.returncode == 0, output[0] + "Bowtie2 build failed\n"
def bowtie2_paired(options):
	b2 = subprocess.Popen(['bowtie2', '-k30', '-t', '-p', options.threads, '-X', options.insert, '-x', options.target, '-1', options.left,  '-2', options.right, '-U', options.unpaired, '-S', options.name], stdout=subprocess.PIPE)	
	output = b2.communicate()
	assert b2.returncode == 0, output[0] + "Bowtie2 alignment failed\n"
def express(options):
	exp = subprocess.Popen(['express', '-o', options.outdir, '-p', options.threads, '%s' % (options.b2base), options.name])	
	output = exp.communicate()
	assert exp.returncode == 0, output[0] + "express failed\n"


##########################################
## alignment depend
##########################################

def checkbowtiebuild():
	try:
	    p = subprocess.Popen(['bowtie-build'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
	    print "Could not find Bowtie2"
	    print "Make sure that it is properly installed on your path"
	    sys.exit(1)
def checkbowtie2():
	try:
	    p = subprocess.Popen(['bowtie2'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
	    print "Could not find Bowtie2"
	    print "Make sure that it is properly installed on your path"
	    sys.exit(1)
def checkexpress():
	try:
	    p = subprocess.Popen(['express'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
	    print "Could not find eXpress"
	    print "Make sure that it is properly installed on your path"
	    sys.exit(1)


##########################################
## Master function
##########################################
def main():
	options = getOptions()
	checkbowtiebuild()	
	checkbowtie2()
	checkexpress()	
	#numbering(options, awker)	
	print >> sys.stderr,"\nBuilding Bowtie2 index, If Necessary: [%s] \n" % (right_now())	
	#b2build(options)
	#print options.target+'.1.bt2'
	if os.path.exists(options.target+'.1.bt2'):
		print >> sys.stderr,"\nLucky You, the Bowtie2 Index Already Exists! I'm going straight to the mapping step. \n"
	else: 
		b2build(options)
	print >> sys.stderr,"\nAligning with Bowtie2: [%s] \n" % (right_now())		
	bowtie2_paired(options)
	print >> sys.stderr,"\nCalculating Expression with eXpress: [%s] \n" % (right_now())
	express(options)
	print >> sys.stderr,"\nDone.. Have a good day! [%s] \n" % (right_now())
if __name__ == "__main__":
	main()

