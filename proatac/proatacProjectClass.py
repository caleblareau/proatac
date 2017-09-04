import click
import os
import os.path
import sys
import shutil
import yaml
import random
import string
import itertools
import time
import platform
from .proatacHelp import *

class proatacProject():
	def __init__(self, script_dir, mode, input, output, name, ncores, bowtie2_index,
		cluster, jobs, keep_duplicates, max_javamem, extract_mito, reference_genome,
		clipl, clipr, keep_temp_files, skip_fastqc,
		bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
		bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path):

		# Figure out operating system
		self.os = "linux"
		if(platform.platform()[0:5]=="Darwi"):
			self.os = "mac"
		
		# verify bowtie2 index
		bwt2idxfiles = os.popen("ls " + bowtie2_index + "*.bt2*").read().strip().split("\n")
		if(len(bwt2idxfiles) < 6):
			sys.exit("ERROR: cannot find bowtie2 index; specify with --bowtie2-index and make sure to add the prefix along with the folder path")
		else:
			self.bowtie2_index = bowtie2_index
		
		# Assign straightforward attributes
		self.clipl = clipl
		self.clipr = clipr
		self.max_javamem = max_javamem
		self.extract_mito = extract_mito
		self.keep_duplicates = keep_duplicates
		self.cluster = cluster
		self.jobs = jobs
		self.name = name
		self.output = output
		self.mode = mode
		
		# Collect samples / fastq lists
		self.samples, self.fastq1, self.fastq2 = inferSampleVectors(input)
		
		self.reference_genome = reference_genome
		# Handle reference genome
		supported_genomes = ['hg19', 'hg38', 'mm9', 'mm10', 'hg19_mm10_c']
		if any(self.reference_genome == s for s in supported_genomes):
			click.echo(gettime() + "Found designated reference genome: %s" % self.reference_genome)
			
			self.tssFile = script_dir + "/anno/TSS/" + self.reference_genome + ".refGene.TSS.bed"
			self.blacklistFile = script_dir + "/anno/blacklist/" + self.reference_genome + ".full.blacklist.bed"
			self.bedtoolsGenomeFile = script_dir + "/anno/bedtools/chrom_" + self.reference_genome + ".sizes"
			
			# Set up effective genome size for macs2
			if self.reference_genome == 'hg19':
				self.BSgenome = 'BSgenome.Hsapiens.UCSC.hg19'
				self.macs2_genome_size = 'hs'
			elif self.reference_genome == 'hg38':
				self.BSgenome = 'BSgenome.Hsapiens.UCSC.hg38'
				self.macs2_genome_size = 'hs'
			elif self.reference_genome == 'mm9':
				self.BSgenome = 'BSgenome.Mmusculus.UCSC.mm9'
				self.macs2_genome_size = 'mm'
			elif self.reference_genome == 'mm10':
				self.BSgenome = 'BSgenome.Mmusculus.UCSC.mm10'
				self.macs2_genome_size = 'mm'
			else:
				self.BSgenome = ''
				self.macs2_genome_size = '4.57e9'
		else: 
			click.echo(gettime() + "Could not identify this reference genome: %s" % self.reference_genome)
			click.echo(gettime() + "Attempting to infer necessary input files from user specification.")
			necessary = [bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome]
			if '' in necessary:
				if reference_genome == '':
					sys.exit("ERROR: specify valid reference genome with --reference-genome flag; QUITTING")
				else:
					sys.exit("ERROR: non-supported reference genome specified so all five of these must be validly specified: --bedtools-genome, --blacklist-file, --tss-file, --macs2-genome-size, --bs-genome; QUITTING")
					
		# Make sure all files are valid
		if(bedtools_genome != ""):
			if(os.path.isfile(bedtools_genome)):
				self.bedtoolsGenomeFile = bedtools_genome
			else: 
				sys.exit("Could not find the bedtools genome file: %s" % bedtools_genome)
				
		if(blacklist_file != ""):
			if(os.path.isfile(blacklist_file)):
				self.blacklistFile = blacklist_file
			else: 
				sys.exit("Could not find the blacklist bed file: %s" % blacklist_file)
		
		if(tss_file != ""):	
			if(os.path.isfile(tss_file)):
				self.tssFile = tss_file
			else: 
				sys.exit("Could not find the bedtools genome file: %s" % tss_file)		
		
		if(bs_genome != ""):
			self.BSgenome = bs_genome
		
		if(macs2_genome_size != ""):
			self.macs2_genome_size = macs2_genome_size
	
		# Setup software paths
		self.bedtools = get_software_path('bedtools', bedtools_path)
		self.bowtie2 = get_software_path('bowtie2', bowtie2_path)
		self.java = get_software_path('java', java_path)
		self.macs2 = get_software_path('macs2', macs2_path)
		self.samtools = get_software_path('samtools', samtools_path)
		self.R = get_software_path('R', r_path)
	