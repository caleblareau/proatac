#!/usr/bin/env python

# Original author: Jason Buenrostro, Stanford University
# Updated by: Caleb Lareau, Harvard University (01 Nov 16)
# The following program will trim PE reads

##### IMPORT MODULES #####
# import necessary for python
import os
import regex
import sys
import bz2
import gzip
import string
import editdistance
import logging
from Bio import SeqIO
from Bio import AlignIO
from optparse import OptionParser

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = editdistance.eval(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-u", action="store_true", default=False, help="Print uncompressed output file")
opts.add_option("-t", default="trim", help="Adds string to trimmed files")
opts.add_option("-p", default="p", help="Prefix for everything")
options, arguments = opts.parse_args()

## Define log but eventually get rid of it
logging.basicConfig(filename='p.log')

# return usage information if no argvs given
if len(sys.argv)==1:
    os.system(sys.argv[0]+" --help")
    sys.exit()

##### INPUTS AND OUTPUTS #####
# name input and outputs
p1_in = options.a
p2_in = options.b

# name outputs and print to working dir
p1_file = options.p + "." + p1_in.split('/')[-1]
p2_file = options.p + "." + p2_in.split('/')[-1]

#check for file type and open input file
append = p1_in.split('.')[-1]
tt = "." + options.t + ".fastq"

if append == "fastq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = regex.sub(".fastq", tt, p1_file)
    p2_out = regex.sub(".fastq", tt, p2_file)
elif append == "fq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = regex.sub(".fq", tt, p1_file)
    p2_out = regex.sub(".fq", tt, p2_file)
elif append == "gz":
    p1_rds = gzip.open(p1_in,'r')
    p2_rds = gzip.open(p2_in,'r')
    p1_out = regex.sub(".fastq.gz", tt, p1_file)
    p2_out = regex.sub(".fastq.gz", tt, p2_file)
elif append == "bz2":
    p1_rds = bz2.BZ2File(p1_in,'r')
    p2_rds = bz2.BZ2File(p2_in,'r')
    p1_out = regex.sub(".fastq.bz2", tt, p1_file)
    p2_out = regex.sub(".fastq.bz2", tt, p2_file)
else:
    sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0;count=1
n=20  # match seq
mismatch=1  # only allow 0-1 mismatches for now

# initilize write files
if options.u == False:
    r1_write = gzip.open(p1_out+'.gz', 'wb')
    r2_write = gzip.open(p2_out+'.gz', 'wb')
elif options.u == True:
    r1_write = open(p1_out, 'w')
    r2_write = open(p2_out, 'w')

while 1:
    # read lines
    p1_line = p1_rds.readline()
    p2_line = p2_rds.readline()

    # break if at end of file
    if not p1_line:
        break

    # load fastq into memory
    if count == 1:
        seqhead1 = p1_line
        seqhead2 = p2_line
    elif count == 2:
        seq1 = p1_line.rstrip()
        seq2 = p2_line.rstrip()
    elif count == 3:
        qualhead1 = p1_line
        qualhead2 = p2_line
    elif count == 4:
        qual1 = p1_line.rstrip()
        qual2 = p2_line.rstrip()

        # align reads to themselves
        i = i+1  # total reads
        rc_seq2 = reverse_complement(seq2[0:n])
        idx = seq1.rfind(rc_seq2) # look for perfect match
        if idx > 0:
            j = j+1  # 0 mismatchs
        elif mismatch>0:
            hold = fuzz_align(rc_seq2,seq1,mismatch)  # else allow for mismatch
            if hold:
                idx,mis=hold
                if mis == 1:
                    k=k+1  # 1 mismatch

        # trim reads if idx exist
        if idx > 0:
            # keep track on how much trimming
            tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
            
            # trim data
            seq1 = seq1[0:idx+n-1] # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
            seq2 = seq2[0:idx+n-1]
            qual1 = qual1[0:idx+n-1]
            qual2 = qual2[0:idx+n-1]
        
        r1_write.write(seqhead1);r1_write.write(seq1+"\n")
        r1_write.write(qualhead1);r1_write.write(qual1+"\n")
        r2_write.write(seqhead2);r2_write.write(seq2+"\n")
        r2_write.write(qualhead2);r2_write.write(qual2+"\n")

    # increment count
    count = count + 1
    if count == 5:
        count = 1
    else:
        count = count

# close files to write the file
r1_write.close();r2_write.close()
p1_rds.close();p2_rds.close()

logging.basicConfig(format='%(asctime)s %(message)s')
logging.info(' Finished trimming data; from the trimmed data--')
logging.info(str(i)+" sequences total")
logging.info(str(j)+" sequences trimmed with 0 mismatches")
logging.info(str(k)+" sequences trimmed with 1 mismatch")
logging.info(str(tot_b/(j+k))+" mean number of bases trimmed for reads requiring trimming. \n")
