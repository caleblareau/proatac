#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University 
# modified from plotV_vC.py (Alicia)

# Will make a V-plot from bed regions

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import re
import pysam
import matplotlib
matplotlib.use('Agg')
from multiprocessing import Pool
from optparse import OptionParser
import random
import scipy.io as sio

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Reads> Accepts sorted BAM file")
opts.add_option("-b", help="<Bed> Accepts bed file")
opts.add_option("-q", default = ".", help="Output directory for log file")
opts.add_option("-s", default = "sample1", help="Sample name")
options, arguments = opts.parse_args()

logoutdir = options.q
sample = options.s

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### DEFINE FUNCTIONS #####
# returns natural sort index
def natural_sort( l ):
	convert = lambda text: int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l,key=alphanum_key) #[headers.index(x) for x in sorted(l,key=alphanum_key)]) ### for index

##### INPUTS AND OUTPUTS #####
# get intervals
p1_ints = np.loadtxt(options.b,delimiter="\t",dtype=bytes).astype(str)

##### SCRIPT #####
# open and read BAM file
bamfile = pysam.Samfile(options.a, "rb")

# initialize vars
outData = {}

# loop through intervals
for i in range(0,len(p1_ints)):
	# initialize rList to avoid double counting
	rList = []
	
	# get interval as num
	s_int=int(p1_ints[i][1])
	e_int=int(p1_ints[i][2])
	
	# loop through rds
	for p2_rds in bamfile.fetch(p1_ints[i][0].tolist(), max(0,s_int-2000), e_int+2000):
		#check mapping quality
		if  p2_rds.qname in rList:
			continue
		
		# get ins pts
		if p2_rds.is_reverse:
			insPos = p2_rds.pos+p2_rds.alen
		else:
			insPos = p2_rds.pos
		
		# if within window
		if insPos > s_int and insPos < e_int:
			# check for barcode
			barcode = dict(p2_rds.tags)['RG']
			if not outData.has_key(barcode):
				outData[barcode] = np.zeros(len(p1_ints))
			
			# add data and save read name
			outData[barcode][i]+=1
			rList.append(p2_rds.qname)
print(outData)
# write output to file
headers = natural_sort(outData.keys())
print(headers)
with open(logoutdir + "/" + sample+ '.tss.log', 'w') as logfile:
	for i in range(0,len(headers)):
		logfile.write(outData[headers[i]])