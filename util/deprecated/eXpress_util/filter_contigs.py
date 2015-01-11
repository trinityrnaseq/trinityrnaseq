#!/usr/bin/env python

#    filter_contigs.py
#    June 2012  Matthew MacManes (macmanes@gmail.com)
#
#   This wrapper is free software: you can redistribute it and/or modify
#


import sys
import subprocess
import optparse
import shutil
import os
from datetime import datetime, date, time
from Bio import SeqIO


print ""
print ""
print ""
print "******************************************************************"
print "*******filter_contigs.py v0.0.1***********************************"
print "*******To run this program, you must have BioPython*******"
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
    parser = optparse.OptionParser(usage="usage: python %prog -b Trinity.fasta  -F [INT]",
                          version="%prog 0.0.1")
    parser.add_option("-b", "--b2base",
                      dest="b2base",
                      default="Trinity.fasta",
                      metavar='file.fa',
                      help="fasta file for B2 index (?Trinity.fasta)")
    parser.add_option("-F", "--FPKM",
                      dest="fpkm",
                      metavar='fpkm',
		      default="1",
                      help="Lowest Acceptable FPKM")
    (options, args) = parser.parse_args()

    return options


file1='results.xprs'
file2='high.exp.results.xprs'
file3='numbers'
##########################################
## alignment procedure
##########################################

def filter_express(options):
	with open('results.xprs','r') as stdin_fh:
		with open('high.exp.results.xprs','w') as stdout_fh:
			num = subprocess.Popen(['awk', '%s >$11{next}1' % (options.fpkm), file1], stdin=stdin_fh, stdout=stdout_fh)
			output = num.communicate()
def make_fasta():
	with open('high.exp.results.xprs','r') as stdin_fh:
		with open('numbers','w') as stdout_fh:
			num = subprocess.Popen(['awk', '{print $2}', file2], stdin=stdin_fh, stdout=stdout_fh)
			output = num.communicate()

def make_fasta_file(options):
	numbers = set()
	with open(file3) as f:
    		for line in f:
        		line = line.strip()
        		if line != "":
        		    try:
        		        num = str(line) #take strings
        		        numbers.add(num)
        		    except:
        		        print line, "Line cannot be parsed"

		fasta_sequences = SeqIO.parse(open('%s' % (options.b2base)),'fasta')
		end = False
		with open('High.cov.%s' % (options.b2base), "w") as f:
	   	 while end != True:
        		try:
        		    seq = fasta_sequences.next()
        		except:
        		    end = True
        		try:
        		    name = str(seq.id) #take strings
        		except:
        		    print seq.id, "Fasta name cannot be parsed"
        		if name in numbers:
        		    SeqIO.write([seq], f, "fasta")
###make fasta






##########################################
## Clean up
##########################################
def clean(options):
	os.remove(file3)

##########################################
## Master function
##########################################
def main():
	options = getOptions()
	filter_express(options)
	make_fasta()	
	make_fasta_file(options)
	clean(options)
	print >> sys.stderr,"\nDone.. Have a good day! [%s] \n" % (right_now())
if __name__ == "__main__":
	main()
