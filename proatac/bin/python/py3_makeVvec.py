#!/usr/bin/env python

# Author: Jason Buenrostro, Stanford University 
# modified from plotV_vC.py (Alicia)

# Will make a V-plot from bed regions

##### IMPORT MODULES #####
# import necessary for python
import os
import sys
import numpy as np
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from multiprocessing import Pool
from optparse import OptionParser
import random

#### OPTIONS ####
# read options from command line
opts = OptionParser()
usage = "usage: %prog [options] [inputs]"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Reads> Accepts sorted BAM file")
opts.add_option("-b", help="<Bed> Accepts bed file")
opts.add_option("-o",help="OutputFile")
opts.add_option("-e",default="1000", help="number of bases to extend to each side, default=1000")
opts.add_option("-p",default="center", help="options:center,ends, default=center")
opts.add_option("-q",help = "Outputplot")
opts.add_option("-c",default="4",help="number of threads to use, default=20")
opts.add_option("-s",default='4',help="column in which strand information is located (1 being first), default=4")
opts.add_option("--window",default='20',help="window size for ploting")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### DEFINE FUNCTIONS #####
# assign mat
def asn_mat(val,mat,s_int,e_int,t,i,weight):
	if float(val)>=s_int and float(val)<e_int-1 and t<rows:  # -1 correct?
		base = val-s_int
		if len(p1_ints[0]) == 3:
			mat[t][base] += weight
		elif p1_ints[i][int(options.s)-1] == "-":
			mat[t][len(mat[0])-base-1] += weight	
		else:
			mat[t][base] += weight
	return mat

#compute Vplot Matrix for a particular chunk
def sub_Mat(start):
	# initialize data matrix
	mat = np.zeros([rows,cols])
	# loop through the intervals and get relevent info
	bamfile = pysam.Samfile(options.a, "rb")
	end=min(start+chunksize,len(p1_ints))
	for i in range(start,end):
		# get interval as num
		center = int(p1_ints[i][1])+(int(p1_ints[i][2])-int(p1_ints[i][1]))//2
		s_int=center-int(options.e)
		e_int=center+int(options.e)
		# loop through rds
		for p2_rds in bamfile.fetch(str(p1_ints[i][0]), max(0,s_int-2000), e_int+2000):
			#check mapping quality
			if p2_rds.mapq<30:# or p2_rds.is_proper_pair==False:
				continue
			# get read positions
			if p2_rds.is_reverse:
				continue
			else:
				l_pos = p2_rds.pos+4
				# calculate center point
				ilen = abs(p2_rds.tlen)-9
				#ilen = 1
				r_pos=l_pos+ilen
				c_pos=l_pos+ilen//2
				if ilen%2==1 and options.p=='center':
					mat=asn_mat(c_pos,mat,s_int,e_int,ilen,i,0.5)
					mat=asn_mat(c_pos+1,mat,s_int,e_int,ilen,i,0.5)
				elif ilen%2!=1 and options.p=='center':
					mat=asn_mat(c_pos,mat,s_int,e_int,ilen,i,1)
					# save ends or read centers to v-plot
				elif options.p == 'ends':
					mat = asn_mat(l_pos,mat,s_int,e_int,ilen,i,1)
					mat = asn_mat(r_pos,mat,s_int,e_int,ilen,i,1)
				else:
					sys.exit('Error, check parameters')
	return mat

##### INPUTS AND OUTPUTS #####
# get intervals
p1_ints = np.loadtxt(options.b, dtype=bytes).astype(str)

##### SCRIPT #####
# open and read BAM file

# determine number of rows and columns for matrix
rows = 1000
cols = int(options.e)*2
#cols = int(p1_ints[0][2])-int(p1_ints[0][1])+int(options.e)*2

# split bedfile into chunks
maxi=len(p1_ints)
chunksize=maxi//int(options.c)
starts=range(0,int(chunksize)*int(options.c)-1,int(chunksize))

# parallel processed computation of matrix for each chunk
if __name__ == "__main__":
	pool = Pool(processes=int(options.c))
	sub_mats=pool.map(sub_Mat, starts, 1)

# sum up matrices for each chunk into matrix for all
mat = np.zeros([rows,cols])
for i in range(len(starts)):
	mat=mat+sub_mats[i]
mat = np.sum(mat,0)

np.savetxt(options.o, np.column_stack((np.arange(-2000,2000),np.array(mat))),delimiter=',',fmt='%s')

fig=plt.figure(figsize=(8.0, 5.0))
xran=min(500,int(options.e))
yran=min(500,rows)
plt.plot(mat/np.mean(mat[1:200]),'k.')
plt.plot(np.convolve(mat,np.ones(int(options.window)),'same')/int(options.window)/np.mean(mat[1:200]),'r')
plt.xlabel('Position to TSS')
plt.ylabel('Insertions')
fig.savefig(options.q)
plt.close(fig)