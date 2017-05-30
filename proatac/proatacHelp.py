import yaml
import itertools
import time

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
	"""
	Basic function to parse yaml/yml
	"""
	samples = []
	if manifest.endswith(('.yaml', '.yml')):
		with open(manifest, 'r') as f: 
			m = yaml.load(f)
		return m
	else:
		click.echo(gettime() + "Please specify a valid .yaml file for analysis")

def gettime(): 
	"""
	Matches `date` in Linux
	"""
	return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
		time.strftime("%Z ") + time.strftime("%Y")+ ": ")
		

def process_seq_dir(d, logf):
	"""
	Function that takes a dictionary parsed from the main .yaml
	and returns something more coherent to be processed downstream
	"""
	name = d['name']
	version = d['version']
	directory = d['dir']
	
	samples = []
	
	if 'remove_samples' in d:
		remove_samples = d['remove_samples']
	else:
		remove_samples = ""
	
	if 'keep_samples' in d:
		keep_samples = d['keep_samples']
	else:
		keep_samples = ""
	
	
	if(keep_samples != ""):
		keeplist = keep_samples.split(",")
	
	if(remove_samples != ""):
		igslist = ignore_samples.split(",")
		for byesample in igslist:
			if byesample in samples: samples.remove(byesample)
	
	
	
	