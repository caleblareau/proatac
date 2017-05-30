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
from .proatacProjectClass import *



# ------------------------------
# Command line arguments
# ------------------------------	

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
	yaml = parse_manifest(manifest)

	
	# -------------------------------
	# Utility functions and variables
	# -------------------------------
	
	outfolder = os.path.abspath(yaml['project_dir']) 
	logfolder = outfolder + "/logs"
	internfolder = outfolder + "/internal"
	parselfolder = internfolder + "/parseltongue"
	
	# Check if directories exist; make if not
	if not os.path.exists(outfolder):
		os.makedirs(outfolder)
	if not os.path.exists(logfolder):
		os.makedirs(logfolder)	
	if not os.path.exists(internfolder):
		os.makedirs(internfolder)
		with open(internfolder + "/README" , 'w') as outfile:
			outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")	
	if not os.path.exists(parselfolder):
		os.makedirs(parselfolder)
		with open(parselfolder + "/README" , 'w') as outfile:
			outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")	
	cwd = os.getcwd()
	
	# Main Project Variable
	p = proatacProject(yaml, script_dir)
	logf = open(logfolder + "/base.proatac.log", 'a')
	
	# -------------------------------
	# Atypical analysis modes
	# -------------------------------	
	
	if check:
		sys.exit("Success! We're reasonably confident that all dependencies and files are good to go!")

	# -----------------------------------
	# ACTUAL PREPROCESSING / SNAKE MAKING
	# -----------------------------------
	
	
	# -------------
	# Adapter Trim
	# -------------
	trimfolder = outfolder + "/01_trimmed"
	if not os.path.exists(trimfolder + "_reads"):
		os.makedirs(trimfolder+ "_reads")
	if not os.path.exists(trimfolder + "_stats"):
		os.makedirs(trimfolder+ "_stats")
				
	click.echo(gettime() + "Trimming samples", logf)
	snakecalla = '''snakemake --snakefile ''' + script_dir + '''/bin/Snakefile.Trim '''
	snakecallb = snakecalla + '''--config allsamples="''' + outfolder + '''/internal/parseltongue/allsamples.csv" '''
	snakecallc = snakecallb + '''outdir="''' + outfolder + '''" scriptdir="''' + script_dir + '''" '''
	snakecall1 = snakecallc + '''python3="''' + p.python3_path + '''" '''
	os.system(snakecall1)
	
	logf.close()
	
	