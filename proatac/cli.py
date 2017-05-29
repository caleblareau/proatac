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
	
		self.yaml = parse_manifest(project_yaml_file_handle)
		self.name = self.yaml['project_name']
		self.project_dir = self.yaml['project_dir']
		self.analysis_person = self.yaml['analysis_person']
		
		# ------------------------------
		# Process reference genome stuff
		# ------------------------------
		
		self.reference_genome = self.yaml['reference_genome']
		
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
	
	project = proatacProject(manifest, script_dir)
	
	# -------------------------------
	# Utility functions and variables
	# -------------------------------
	
	outfolder = os.path.abspath(project.project_dir) 
	logfolder = outfolder + "/logs"
	 
	cwd = os.getcwd()
	def gettime(): # Matches `date` in Linux
		return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
			time.strftime("%Z ") + time.strftime("%Y")+ ": ")
	
	# -------------------------------
	# Atypical analysis modes
	# -------------------------------	
	
	if check:
		sys.exit("ERROR: Haven't actually set this up yet")

	# -------------------------------
	# Now actually do stuff
	# -------------------------------

	# Check if directory exists; make it if not
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	
	
	
	print(gettime())
	
	
	