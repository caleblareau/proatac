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
from pkg_resources import get_distribution
from subprocess import call, check_call
from .proatacHelp import *

# -----------------------
# Core project option
# -----------------------

class proatacProject():
	def __init__(self, project_yaml_file_handle, script_dir):
		
		# Basic attributes
		self.yaml = parse_manifest(project_yaml_file_handle)
		self.name = self.yaml['project_name']
		self.project_dir = self.yaml['project_dir']
		self.analysis_person = self.yaml['analysis_person']
		
		# ------------------------------
		# Make folders / log files
		# ------------------------------
		
		outfolder = os.path.abspath(self.project_dir) 
		logfolder = outfolder + "/logs"
		
		# Check if directories exist; make if not
		if not os.path.exists(outfolder):
			os.makedirs(outfolder)
		if not os.path.exists(logfolder):
			os.makedirs(logfolder)
		logf = open(logfolder + "/base.proatac.log", 'a')
		
		# ------------------------------
		# Process reference genome stuff
		# ------------------------------
		
		self.reference_genome = self.yaml['reference_genome']
		supported_genomes = ['hg19', 'hg38', 'mm9', 'mm10', 'hg19_mm10']
		if any(self.reference_genome in s for s in supported_genomes):
			click.echo(gettime() + "Found designated reference genome: %s" % self.reference_genome, logf)
		else: 
			click.echo(gettime() + "Could not identify this reference genome: %s" % self.reference_genome, logf)
		
		# ------------------------
		# Process dependency paths
		# ------------------------
		
		# bowtie2
		if(self.yaml['paths']['bowtie2_path'] != ''):
			self.bowtie2_path = self.yaml['paths']['bowtie2_path']
		else:
			self.bowtie2_path = shutil.which("bowtie2")
		if(self.bowtie2_path == "None"):
			sys.exit("ERROR: cannot find bowtie2 in environment; set the 'bowtie2_path' in the .yaml file or add to PATH")
		
		# macs2	
		if(self.yaml['paths']['macs2_path'] != ''):
			self.macs2_path = self.yaml['paths']['macs2_path']
		else:
			self.macs2_path = shutil.which("macs2")
		if(self.macs2_path == "None"):
			sys.exit("ERROR: cannot find macs2 in environment; set the 'macs2_path' in the .yaml file or add to PATH")
		
		# samtools
		if(self.yaml['paths']['samtools_path'] != ''):
			self.samtools_path = self.yaml['paths']['samtools_path']
		else:
			self.samtools_path = shutil.which("samtools")
		if(self.samtools_path == "None"):
			sys.exit("ERROR: cannot find samtools in environment; set the 'samtools_path' in the .yaml file")
				
		# R
		if(self.yaml['paths']['R_path'] != ''):
			self.R_path = self.yaml['paths']['R_path']
		else:
			self.R_path = shutil.which("R")
		if(self.R_path == "None"):
			sys.exit("ERROR: cannot find R in environment; set the 'R_path' in the .yaml file")
						
		# Python3
		if(self.yaml['paths']['python3_path'] != ''):
			self.python3_path = self.yaml['paths']['python3_path']
		else:
			self.python3_path = shutil.which("python3")
		if(self.python3_path == "None"):
			sys.exit("ERROR: cannot find python3 in environment; set the 'python3_path' in the .yaml file")
		
		# Check for R package dependencies
		required_packages = ['ggplot2', 'tidyverse']
		installed_packages = os.popen('''R -e "installed.packages()" | awk '{print $1}' | sort | uniq''').read().strip().split("\n")
		if(not set(required_packages) < set(installed_packages)):
			sys.exit("ERROR: cannot find the following R package: " + str(set(required_packages) - set(installed_packages)) + "\n" + 
				"Install it in your R console and then try rerunning proatac.")
		
		# ------------------------------
		# Process sequencing directories
		# ------------------------------
				
		for run in self.yaml['sequencing_directories']:
			print(run)
		
		

@click.command()
@click.version_option()
@click.option('--check', is_flag=True, help='[MODE] Check to see if all dependencies are properly configured')
@click.option('--stingy', is_flag=True, help='Space-efficient analyses; remove ')
@click.argument('manifest')

def main(manifest, check, stingy):
	"""Preprocessing ATAC and scATAC Data."""
	__version__ = get_distribution('proatac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting proatac pipeline v%s" % __version__)
	project = proatacProject(manifest, script_dir)
	
	# -------------------------------
	# Utility functions and variables
	# -------------------------------
	
	outfolder = os.path.abspath(project.project_dir) 
	logfolder = outfolder + "/logs"
	logf = open(logfolder + "/base.proatac.log", 'a')
	cwd = os.getcwd()
	
	# -------------------------------
	# Atypical analysis modes
	# -------------------------------	
	
	if check:
		sys.exit("ERROR: Haven't actually set this up yet")

	# -------------------------------
	# Now actually do stuff
	# -------------------------------



	
	logf.close()
	
	