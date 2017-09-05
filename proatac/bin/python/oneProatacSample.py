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
name = config["name"]
script_dir = config["script_dir"]
mode = config["mode"]

# Ned to update
keepchrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
mitochrs = ['chrM', 'MT', 'humanM', 'mouseM', 'humanMT', 'mouseMT']
read_quality = "20"
# bigwig support

# Step 1
clipl = config["clipl"]
clipr = config["clipr"]
skip_fastqc = str(config["skip_fastqc"])
java = str(config["java"])

# Step 2
bowtie2 = config["bowtie2"]
samtools = config["samtools"]
bowtie2index = config["bowtie2_index"]
extract_mito = str(config["extract_mito"])

# Step 3
max_javamem = config["max_javamem"]
keep_duplicates = str(config["keep_duplicates"])
MarkDuplicatesCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/MarkDuplicates.jar"

# Step 4
macs2 = config["macs2"]
macs2_genome_size = config["macs2_genome_size"]

# 01 Trim using custom script
trim_py = script_dir + "/bin/python/py3_ATACtrim.py"
pycall = "python " + trim_py + " -a "+fastq1+" -b "+fastq2+" -l "+clipl+" -r "+clipr+" -s "+sample+" -o "+outdir+"/01_trimmed -t hard -q "+outdir+"/logs/trim"
if not os.path.isfile(outdir+"/logs/trim/"+sample+".trim.log"):
	os.system(pycall)
tfq1 = outdir + "/01_trimmed/" + sample + "_1.trim.fastq.gz"
tfq2 = outdir + "/01_trimmed/" + sample + "_2.trim.fastq.gz"

# 01a fastqc
fastqc_path = shutil.which("fastqc")
if(skip_fastqc == "False" and str(fastqc_path) != "None"):
	fastqc_call = "fastqc -q -j "+ java + " " + tfq1 + " " + tfq2 + " -o "+outdir+"/logs/fastqc" 
	os.system(fastqc_call)
	
# 02 Align with bowtie2
bwt2log = outdir+"/logs/bowtie2/" + sample + ".log"
sortedbam = outdir + '/02_aligned_reads/'+sample+'.all.sorted.bam'
bwt2call = '(' + bowtie2 + ' -X 2000 -p 2 -x '+ bowtie2index +' --rg-id '+sample+' -1 '+tfq1+' -2 '+tfq2+' | ' + samtools + ' view -bS - | ' + samtools + ' sort -@ 2 - -o ' +sortedbam+') 2> ' + bwt2log 
if not os.path.isfile(sortedbam):
	os.system(bwt2call)
	pysam.index(sortedbam)

# 02a mito-- extract if the user says to
if(extract_mito == "True"):
	mitobam = outdir + "/mito/" + sample + ".mito.bam"
	mitocall = samtools + " view -b "+sortedbam+" -o "+mitobam+" " + " ".join(str(i) for i in mitochrs) + " && " + samtools + " index " + mitobam
	if not os.path.isfile(mitobam):
		os.system(mitocall)

# 03 Process .bam files
temp1 = outdir + "/03_processed_reads/temp/" + sample + ".temp1.bam"
temp1call = samtools + " view -q " + read_quality + " -f 0x2 "+sortedbam+" -o "+ temp1
if not os.path.isfile(temp1):
	os.system(temp1call)
	pysam.index(temp1)

temp2 = outdir + "/03_processed_reads/temp/" + sample + ".temp2.bam"
temp2call = samtools + " view -b "+temp1+" -o "+temp2+" " + " ".join(str(i) for i in keepchrs)
if not os.path.isfile(temp2):
	os.system(temp2call)
	pysam.index(temp2)

if(mode == "bulk"):
	finalbam = outdir + "/final/bams/"+sample+".proatac.bam"
else:
	finalbam = outdir + "/03_processed_reads/bams/"+sample+".proatac.bam"	
rmlog = outdir + "/logs/picard/markdups/"+sample+".MarkDuplicates.log"
run_markdup = MarkDuplicatesCall + " INPUT="+temp2+" OUTPUT="+finalbam+" METRICS_FILE="+rmlog+" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT"

if not os.path.isfile(finalbam):
	if(keep_duplicates == "True"):
		os.system("cp " + temp2 + " " + finalbam)
	else:
		os.system(run_markdup)
	pysam.index(finalbam)


# 04 Do macs2 on each sample only if bulk
# Need to build in support for bigwig here
if(mode == "bulk"):
	macs2log = outdir + "/logs/macs2/"+sample+".peakcalls.log"
	macs2outdir = outdir + "/04_qc/macs2_each/"
	macs2call = "(" + macs2 + " callpeak -t "+finalbam+" -n " + macs2outdir + sample + " --nomodel --keep-dup all --call-summits -q 0.05 -g " + macs2_genome_size + ") 2> " + macs2log
	if(macs2outdir + sample + "_peaks.narrowPeak"):
		os.system(macs2call)
		os.system("mv " + macs2outdir + sample + "_peaks.xls " + outdir + "/logs/macs2/"+sample+"_peaks.xls")
		os.system("mv " + macs2outdir + sample + "_summits.bed " + outdir + "/final/summits/"+sample+"_summits.bed")

os.system('echo "hey" > ' + outdir+"/logs/trim/"+sample+".trim.txt")