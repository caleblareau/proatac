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

# ----------------------------------
# Core object for handling the .yaml
# ----------------------------------

class proatacProject():
	def __init__(self, yaml, script_dir):
		
		# Basic attributes
		self.yaml = yaml
		self.project_name = self.yaml['project_name']
		self.project_dir = self.yaml['project_dir']
		self.analysis_person = self.yaml['analysis_person']
		
		if((' ' in self.yaml['project_name']) or (' ' in self.yaml['project_dir'])):
			sys.exit("ERROR: remove white space from project_name and project_dir variables in .yaml")
		
		# ------------------------------
		# Process reference genome stuff
		# ------------------------------
		
		# Computing configuration
		if ("read_quality" in self.yaml['parameters']) and (str(self.yaml['parameters']['read_quality']) != "None"):
			if(int(self.yaml['parameters']['read_quality']) < 0):
				sys.exit("ERROR: read_quality as specified in the .yaml seems fishy; make sure it's a positive integer")
			self.read_quality = str(int(self.yaml['parameters']['read_quality']))
		else:
			self.read_quality = str(30)

		# Computing configuration
		if ("max_cores" in self.yaml['parameters']) and (str(self.yaml['parameters']['max_cores']) != "None"):
			self.max_cores = str(self.yaml['parameters']['max_cores'])
		else:
			self.max_cores = str(1)
		
		# Computing configuration
		if ("extract_mito" in self.yaml['parameters']) and (str(self.yaml['parameters']['extract_mito']) != "None"):
			self.extract_mito = str(self.yaml['parameters']['extract_mito'])
		else:
			self.extract_mito = "false"

		# Computing configuration
		if ("chr_name_length" in self.yaml['parameters']) and (str(self.yaml['parameters']['chr_name_length']) != "None"):
			self.chr_name_length = str(self.yaml['parameters']['chr_name_length'])
		else:
			self.chr_name_length = str(10)
		
		# Figure out operating system
		self.os = "linux"
		if(platform.platform()[0:5]=="Darwi"):
			self.os = "mac"
		
		if(self.os == "mac"):
			self.peat_path = script_dir + "/bin/mac/PEAT_mac"
			self.pigz_path = script_dir + "/bin/mac/pigz_mac"
		else:
			self.peat_path = script_dir + "/bin/linux/PEAT_linux"
			self.pigz_path = script_dir + "/bin/linux/pigz_linux"
		
		outfolder = os.path.abspath(yaml['project_dir']) 
		logfolder = outfolder + "/logs"
		logf = open(logfolder + "/base.proatac.log", 'a')
		
		# ------------------------------
		# Process peak information
		# ------------------------------
		
		self.reference_genome = self.yaml['reference_genome']
		supported_genomes = ['hg19', 'hg38', 'mm9', 'mm10', 'hg19_mm10_c']
		if any(self.reference_genome in s for s in supported_genomes):
			click.echo(gettime() + "Found designated reference genome: %s" % self.reference_genome, logf)
			self.tssFile = script_dir + "/anno/TSS/" + self.reference_genome + ".refGene.TSS.bed"
			self.blacklistFile = script_dir + "/anno/blacklist/" + self.reference_genome + ".ful.blacklist.bed"
			
			# Set up effective genome size for macs2
			if any(self.reference_genome in s for s in ['hg19', 'hg38']):
				self.macs2_genome_size = 'hs'
			else if any(self.reference_genome in s for s in ['mm9', 'mm10']):
				self.macs2_genome_size = 'mm'
			else:
				self.macs2_genome_size = '4.57e9'
		else: 
			click.echo(gettime() + "Could not identify this reference genome: %s" % self.reference_genome, logf)
				
		if ("peak_settings" in self.yaml['parameters']):
			if "tssFile" in self.yaml['peak_settings']:
				b = self.yaml['peak_settings']['tssFile']
				if(b != ''):
					self.tssFile = os.path.realpath(b)
			if "blacklistFile" in self.yaml['peak_settings']:
				b = self.yaml['peak_settings']['blacklistFile']
				if(b != ''):
					self.blacklistFile = os.path.realpath(b)


		# ------------------------
		# Process dependency paths
		# ------------------------
		
		# bowtie2
		if(self.yaml['paths']['bowtie2_path'] != ''):
			self.bowtie2_path = self.yaml['paths']['bowtie2_path']
		else:
			self.bowtie2_path = shutil.which("bowtie2")
		if(str(self.bowtie2_path) == "None"):
			sys.exit("ERROR: cannot find bowtie2 in environment; set the 'bowtie2_path' in the .yaml file or add to PATH")
		
		# bowtie2 index
		bwt2idxfiles = os.popen("ls " + self.yaml['paths']['bowtie2_index']+ "*.bt2").read().strip().split("\n")
		if(len(bwt2idxfiles) < 6):
			sys.exit("ERROR: cannot find bowtie2 index; make sure to add the prefix along with the folder path")
		else:
			self.bowtie2_index = self.yaml['paths']['bowtie2_index']
		
		# macs2	
		if(self.yaml['paths']['macs2_path'] != ''):
			self.macs2_path = self.yaml['paths']['macs2_path']
		else:
			self.macs2_path = shutil.which("macs2")
		if(str(self.macs2_path) == "None"):
			sys.exit("ERROR: cannot find macs2 in environment; set the 'macs2_path' in the .yaml file or add to PATH")
		
		# samtools
		if(self.yaml['paths']['samtools_path'] != ''):
			self.samtools_path = self.yaml['paths']['samtools_path']
		else:
			self.samtools_path = shutil.which("samtools")
		if(str(self.samtools_path) == "None"):
			sys.exit("ERROR: cannot find samtools in environment; set the 'samtools_path' in the .yaml file or add to PATH")
				
		# R
		if(self.yaml['paths']['R_path'] != ''):
			self.R_path = self.yaml['paths']['R_path']
		else:
			self.R_path = shutil.which("R")
		if(str(self.R_path) == "None"):
			sys.exit("ERROR: cannot find R in environment; set the 'R_path' in the .yaml file or add to PATH")

		# Java
		if(self.yaml['paths']['java_path'] != ''):
			self.java_path = self.yaml['paths']['java_path']
		else:
			self.java_path = shutil.which("java")
		if(str(self.java_path) == "None"):
			sys.exit("ERROR: cannot find Java in environment; set the 'java_path' in the .yaml file or add to PATH")
			
		# bedtools
		if(self.yaml['paths']['bedtools_path'] != ''):
			self.java_path = self.yaml['paths']['bedtools_path']
		else:
			self.bedtools_path = shutil.which("bedtools")
		if(str(self.bedtools_path) == "None"):
			sys.exit("ERROR: cannot find bedtools in environment; set the 'bedtools_path' in the .yaml file or add to PATH")
												
		# Python3
		if(self.yaml['paths']['python3_path'] != ''):
			self.python3_path = self.yaml['paths']['python3_path']
		else:
			self.python3_path = shutil.which("python3")
		if(str(self.python3_path) == "None"):
			sys.exit("ERROR: cannot find python3 in environment; set the 'python3_path' in the .yaml file or add to PATH")
		
		# Check for R package dependencies
		required_packages = ['ggplot2', 'tidyverse']
		installed_packages = os.popen(self.R_path + ''' -e "installed.packages()" | awk '{print $1}' | sort | uniq''').read().strip().split("\n")
		if(not set(required_packages) < set(installed_packages)):
			sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
				"Install it in your R console and then try rerunning proatac (but there may be other missing dependencies).")
		
		# The final step should be fast, so remove the file that coordinates all samples if it exists
		listAllSamples = outfolder + '/.internal/parseltongue/allsamples.csv'
		if os.path.exists(listAllSamples):
			os.remove(listAllSamples)
		
		# Process sequencing directories		
		for run in self.yaml['sequencing_directories']:
			process_seq_dir(run, logf, listAllSamples)
		
		# Check to see if any sample names are duplicated
		with open(listAllSamples) as f:
			seen = set()
			for line in f:
				line_lower = line.split(",")[0]
				if line_lower in seen:
					sys.exit("ERROR: found multiple sample IDs specified the same way; this will cause problems down the road; quitting now")
				else:
					seen.add(line_lower)
		
		logf.close()
