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
			else:
				# Bedtools genome file
				if(os.path.isfile(bedtools_genome)):
					self.bedtoolsGenomeFile = bedtools_genome
				else: 
					sys.exit("Could not find the bedtools genome file: %s" % bedtools_genome)

				# Blacklist coordinates file
				if(os.path.isfile(blacklist_file)):
					self.blacklistFile = blacklist_file
				else: 
					sys.exit("Could not find the blacklist bed file: %s" % blacklist_file)
					
				# TSS coordinates file
				if(os.path.isfile(tss_file)):
					self.tssFile = tss_file
				else: 
					sys.exit("Could not find the bedtools genome file: %s" % tss_file)		
				self.BSgenome = bs_genome
				self.macs2_genome_size = macs2_genome_size

