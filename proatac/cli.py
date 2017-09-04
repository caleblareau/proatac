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
import pysam

from pkg_resources import get_distribution
from subprocess import call, check_call
from .proatacHelp import *
from .proatacProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['bulk', 'check', 'single', 'summitsToPeaks']))

@click.option('--input', '-i', default = ".", required=True, help='Input for particular mode; see documentation')
@click.option('--output', '-o', default="proatac_out", help='Output directory for analysis; see documentation.')
@click.option('--name', '-n', default="proatac",  help='Prefix for project name')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--bowtie2-index', '-bi', default = "", required=True, help='Path to the bowtie2 index; should be specified as if you were calling bowtie2 (with file index prefix)')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--peak-width', '-pw', default = "250", help='Fixed width value of resulting peaks from summit calling / padding. 250, 500 recommended.')
@click.option('--keep-duplicates', '-kd', is_flag=True, help='Keep optical/PCR duplicates')
@click.option('--max-javamem', '-jm', default = "4000m", help='Maximum memory for java')
@click.option('--extract-mito', '-em', is_flag=True, help='Extract mitochondrial reads as part of special output.')
@click.option('--reference-genome', '-rg', default = "", help='Support for built-in genome; choices are hg19, mm9, hg38, mm10, hg19_mm10_c (species mix)')

@click.option('--clipL', '-cl', default = "0", help='Number of variants to clip from left hand side of read.')
@click.option('--clipR', '-cr', default = "0", help='Number of variants to clip from right hand side of read.')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--skip-fastqc', '-sf', is_flag=True, help='Throw this flag to skip fastqc on the trimmed fastqs; will only run if software is discovered in the path')

@click.option('--bedtools-genome', '-bg', default = "", help='Path to bedtools genome; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--blacklist-file', '-bl', default = "", help='Path to bed file of blacklist; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--tss-file', '-ts', default = "", help='Path bed file of transcription start sites; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--macs2_genome_size', '-mg', default = "", help='String passed to macs2 for ; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--bs-genome', '-bs', default = "", help='String corresponding to the R/Bioconductor package for BS genome of build; overrides default if --reference-genome flag is set and is necessary for non-supported genomes..')

@click.option('--bedtools-path', default = "", help='Path to bedtools; by default, assumes that bedtools is in PATH')
@click.option('--bowtie2-path', default = "", help='Path to bowtie2; by default, assumes that bowtie2 is in PATH')
@click.option('--java-path', default = "", help='Path to java; by default, assumes that java is in PATH')
@click.option('--macs2-path', default = "", help='Path to macs2; by default, assumes that macs2 is in PATH')
@click.option('--samtools-path', default = "", help='Path to samtools; by default, assumes that samtools is in PATH')
@click.option('--R-path', default = "", help='Path to R; by default, assumes that R is in PATH')

def main(mode, input, output, name, ncores, bowtie2_index,
	cluster, jobs, peak_width, keep_duplicates, max_javamem, extract_mito, reference_genome,
	clipl, clipr, keep_temp_files, skip_fastqc,
	bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
	bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path):
	
	"""
	proatac: a toolkit for prePROcessing ATAC-seq and scatac-seq data. \n
	Caleb Lareau, Aryee/Buenrostro Labs \n
	modes = ['bulk', 'check', 'single'] \n
	See http://proatac.readthedocs.io for more details.
	"""
	
	__version__ = get_distribution('proatac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting proatac pipeline v%s" % __version__)
	
	# Take a collection of summits files and return a consensus set of peaks
	# Based on the relative ranking approach	
	if(mode == 'summitsToPeaks'):
		make_folder(
		doSummitsToPeaks(input, name, peak_width)
		click.echo(gettime() + "Starting proatac pipeline v%s"
		
	# Make a mode to handle split-pool data
	
	p = proatacProject(script_dir, mode, input, output, name, ncores, bowtie2_index,
		cluster, jobs, peak_width, keep_duplicates, max_javamem, extract_mito, reference_genome,
		clipl, clipr, keep_temp_files, skip_fastqc,
		bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
		bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path)
	
	if (mode == "check"):
		click.echo(gettime() + "Dependencies and user-reported file paths OK")
		click.echo("\nproatac will process the following samples / files with bulk / single specified: \n")
		print("Sample", "Fastq1", "Fastq2")
		for x in range(len(p.samples)):
			print(p.samples[x], p.fastq1[x], p.fastq2[x])
		click.echo("\nIf this table doesn't look right, consider specifying a manually created sample input table (see documentation).\n")
		sys.exit(gettime() + "Successful check complete; QUITTING.")
	
	# Single or bulk processing
	if(mode == "single" or mode == "bulk"):
	
		of = output; logs = of + "/logs"; fin = of + "/final"; trim = of + "/01_trimmed"; 
		aligned = of + "/02_aligned_reads"; processed = of + "/03_processed_reads";
		qc = of + "/04_qc"
		
		folders = [of, logs, fin, trim, aligned, processed, qc,
			of + "/.internal/parseltongue", of + "/.internal/samples"]

		mkfolderout = [make_folder(x) for x in folders]
		
		# Determine chromosomes to keep / filter
		#chrs = os.popen(p.samtools_path + " idxstats " +  outfolder + "/02_aligned_reads/* | cut -f1").read().strip().split("\n")
		rmchrlist = ["*", "chrY", "MT", "chrM"]
		#keepchrs = [x for x in chrs if x not in rmchrlist and len(x) < int(p.chr_name_length)]
		
		# To do: handle snakemake file a la mgatk -- one scatter (probably) file; one gather?
		# Create a sample table / dump it to internal folder
		# Dump the project to a .yaml file for passing around (a la mgatk)
		
		if keep_temp_files:
			click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.")
		else:
			if(mode == "bulk" or mode == "single"):
				byefolder = of
			
			shutil.rmtree(byefolder + "/logs")
			shutil.rmtree(byefolder + "/.internal")
			shutil.rmtree(byefolder + "/01_trimmed")
			shutil.rmtree(byefolder + "/02_aligned_reads")
			shutil.rmtree(byefolder + "/03_processed_reads")
			shutil.rmtree(byefolder + "/04_qc")
			click.echo(gettime() + "Intermediate files successfully removed.")
		
	click.echo(gettime() + "Complete.")