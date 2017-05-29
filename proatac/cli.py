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
        self.analyst = self.yaml['analysis_person']
		
        self.libraries = OrderedDict()
        self.runs = OrderedDict()

        for run in self.yaml['sequencing_runs']:
            print(run)

    @property
    def parameters(self):
        if not hasattr(self, '_parameters'):
            #Read defaults
            with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default.yaml'), 'r') as f:
                self._parameters = yaml.load(f)['parameters']
            # Update with user provided values
            if 'parameters' in self.yaml:
                for k, d in self.yaml['parameters'].items():
                    self._parameters[k].update(d)

        return self._parameters


    @property
    def paths(self):
        if not hasattr(self, '_paths'):
            script_dir = os.path.dirname(os.path.realpath(__file__))
            #Read defaults
            with open(os.path.join(script_dir, 'default_parameters.yaml'), 'r') as f:
                paths = yaml.load(f)['paths']
            # Update with user provided values
            paths.update(self.yaml['paths'])

            paths['python'] = os.path.join(paths['python_dir'], 'python')
            paths['macs2'] = os.path.join(paths['macs2_dir'], 'macs2')
            paths['bowtie'] = os.path.join(paths['bowtie_dir'], 'bowtie')
            paths['samtools'] = os.path.join(paths['samtools_dir'], 'samtools')
            paths['rsem_tbam2gbam'] = os.path.join(paths['rsem_dir'], 'rsem-tbam2gbam')
            paths['rsem_prepare_reference'] = os.path.join(paths['rsem_dir'], 'rsem-prepare-reference')


@click.command()
@click.option('--check', default="hichipper_out", required=False, help='Throw this flag to check if all dependencies are configured in the .yaml')
@click.argument('manifest')
@click.version_option()

def main(manifest, check):
	"""Preprocessing ATAC and scATAC Data."""
	__version__ = get_distribution('proatac').version
	
	def gettime(): # Matches `date` in Linux
		return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + time.strftime("%Z ") + time.strftime("%Y")+ ": ")
	
	if check:
		sys.exit("ERROR: Output path (%s) already exists." % out)
	

	print(gettime())
	project = proatacProject(manifest)
	
	