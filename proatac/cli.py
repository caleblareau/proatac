import click
import os
import os.path
import sys
import shutil
import yaml
import shutil
import random
import string
import itertools
import time
from pkg_resources import get_distribution
from subprocess import call, check_call



# -----------------------
# Helper functions
# -----------------------

def string_hamming_distance(str1, str2):
    """
    Fast hamming distance over 2 strings known to be of same length.
    In information theory, the Hamming distance between two strings of equal 
    length is the number of positions at which the corresponding symbols 
    are different.
    eg "karolin" and "kathrin" is 3.
    """
    return sum(itertools.imap(operator.ne, str1, str2))

def rev_comp(seq):
    """
    Fast Reverse Compliment
    """  
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def parse_manifest(manifest):
    samples = []
    if manifest.endswith(('.yaml', '.yml')):
        with open(manifest, 'r') as f: 
            m = yaml.load(f)
        return m
    else:
        click.echo(gettime() + "Please specify a valid .yaml file for samples")


# -----------------------
# Core objects
# -----------------------

class proatacProject():
	def __init__(self, project_yaml_file_handle):
		self.yaml = parse_manifest(project_yaml_file_handle)
		print(self.yaml['paths'])
		self.name = self.yaml['project_name']
		self.project_dir = self.yaml['project_dir']
		self.analysis_person = self.yaml['analysis_person']
		self.reference_genome = self.yaml['reference_genome']
		
		for run in self.yaml['sequencing_directories']:
			print(run)


@click.command()
@click.option('--check', is_flag=True, help='[MODE] Check to see if all dependencies are properly configured')
@click.argument('manifest')
@click.version_option()

def main(manifest, check):
	"""Preprocessing ATAC and scATAC Data."""
	__version__ = get_distribution('proatac').version
	
	def gettime(): # Matches `date` in Linux
		return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + time.strftime("%Z ") + time.strftime("%Y")+ ": ")
	
	project = proatacProject(manifest)
	
	if check:
		sys.exit("ERROR: Haven't actually set this up yet")
	
	# Check if directory exists; make it if not
	if not os.path.exists(project.project_dir):
		os.makedirs(project.project_dir)

	print(gettime())
	
	
	