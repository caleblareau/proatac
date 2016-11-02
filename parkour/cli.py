import click
import os
import os.path
import sys
import shutil
import shutil
import random
import string
import logging
from pkg_resources import get_distribution
from subprocess import call, check_call

#Grab Sample Names
def get_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

@click.command()
@click.argument('mode')
@click.option('-o', default="p", help='Output prefix associated with sample')
@click.option('-a', default="", help='Filename/path for read 1 for sample')
@click.option('-b', default="", help='Filename/path for read 2 for sample')
@click.option('-u', is_flag=True, help='Leave output .fastq files uncompressed?')
@click.option('-s', default="-p 0.01 --nomodel", help='String of arguments to pass to MACS2; default = "-p 0.01 --nomodel"')
@click.option('-q', is_flag=True, help='Skip QC report generation? (Requires R + dependent packages (see README))')
@click.version_option()



def main(mode, o, a, b, u, s, q):
	
##################
#### MODES I/O ###
##################

	"""
	Valid MODE options include `trim`,
	
	\b
	`trim` mode valid options:
	  -a file 1
	  -b file 2
	
	
	"""
	__version__ = get_distribution('parkour').version
	modes = ['trim']
	if not any(mode in s for s in modes):
		sys.exit("ERROR: Improper mode '" + mode + "' selected")
	
	click.echo("Running " + mode +" mode in parkour v%s" % __version__)
	script_dir = os.path.dirname(os.path.realpath(__file__))
	

	
###################
#### TRIM MODE ####
###################

	if mode == 'trim':
		# Check that all parameters that are needed are here
		if a == "":
			sys.exit("ERROR: Supply an argument with -a to run trim mode")
		if b == "":
			sys.exit("ERROR: Supply an argument with -b to run trim mode")
		
		uncmprs = str(False)
		if u:
			uncmprs = str(True)
		
		# Check that all the parameters are valid
		if not os.path.isfile(a):
			sys.exit("ERROR: File '" + a + "' specified with -a does not exist!")
		if not os.path.isfile(b):
			sys.exit("ERROR: File '" + b + "' specified with -b does not exist!")

		cmd = ['python', os.path.join(script_dir, 'pyadapter_trim.py'), "-a", str(a), "-b", str(b), "-u", str(uncmprs), "-o", str(o)] 
		click.echo(cmd)
		call(cmd)
	click.echo("Done")