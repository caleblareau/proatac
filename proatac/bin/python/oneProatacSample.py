#!/usr/bin/env python

from os.path import join
import os
import subprocess
import sys
import shutil
import pysam
from ruamel import yaml

configFile = sys.argv[1]
fastq1 = sys.argv[2]
fastq2 = sys.argv[3]
sample = sys.argv[4]

with open(configFile, 'r') as stream:
	config = yaml.load(stream, Loader=yaml.Loader)

# General
outdir = config["output"]
script_dir = config["script_dir"]

# Step 1
clipl = config["clipl"]
clipr = config["clipr"]
skip_fastqc = str(config["skip_fastqc"])

# Step 2
bowtie2_index = config["bowtie2_index"]

# 01 Trim
trim_py = script_dir + "/bin/python/py3_ATACtrim.py"
pycall = "python " + trim_py + " -a "+fastq1+" -b "+fastq2+" -l "+clipl+" -r "+clipr+" -s "+sample+" -o "+outdir+"/01_trimmed -t hard -q "+outdir+"/logs/trim"
os.system(pycall)

tfq1 = outdir + "/01_trimmed/" + sample + "_1.trim.fastq.gz"
tfq2 = outdir + "/01_trimmed/" + sample + "_2.trim.fastq.gz"

# 01a fastqc
fastqc_path = shutil.which("fastqc")
if(skip_fastqc == "false" and str(fastqc_path) != "None"):
	fastqc_call = "fastqc " + tfq1 + " " + tfq2 + " -o "+outdir+"/logs/fastqc" 
	print(fastqc_call)
	os.system(fastqc_call)
	
# 02 Align
#bwt2call = '(' + bowtie2 + ' -X 2000 -p {threads} -x '+ bowtie2index +' --rg-id {wildcards.sample} -1 {input.r1} -2 {input.r2} | ' + samtools + ' view -bS - | ' + samtools + ' sort -@ {threads} - -o ' + outdir + '/02_aligned_reads/{wildcards.sample}.all.sorted.bam) 2> {log}' 
