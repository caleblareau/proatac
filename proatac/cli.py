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

@click.argument('mode', type=click.Choice(['bulk', 'check', 'single']))

@click.option('--input', '-i', default = ".", required=True, help='Input for particular mode; see documentation')
@click.option('--output', '-o', default="proatac_out", help='Output directory for analysis; see documentation.')
@click.option('--name', '-n', default="proatac",  help='Prefix for project name')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--bowtie2-index', '-bi', default = ".", required=True, help='Path to the bowtie2 index; should be specified as if you were calling bowtie2 (with file index prefix)')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

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
@click.option('--tss-file', '-ts', default = "", help='Path bed file of transcription start sites; overrides default if --reference-genome flag is set and is necessary for non-supported genomes..')
@click.option('--macs2_genome_size', '-mg', default = "", help='String passed to macs2 for ; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--bs-genome', '-bs', default = "", help='String corresponding to the R/Bioconductor package for BS genome of build; overrides default if --reference-genome flag is set and is necessary for non-supported genomes..')

@click.option('--bedtools-path', default = "", help='Path to bedtools; by default, assumes that bedtools is in PATH')
@click.option('--bowtie2-path', default = "", help='Path to bowtie2; by default, assumes that bowtie2 is in PATH')
@click.option('--java-path', default = "", help='Path to java; by default, assumes that java is in PATH')
@click.option('--macs2-path', default = "", help='Path to macs2; by default, assumes that macs2 is in PATH')
@click.option('--samtools-path', default = "", help='Path to samtools; by default, assumes that samtools is in PATH')
@click.option('--R-path', default = "", help='Path to R; by default, assumes that R is in PATH')

def main(mode, input, output, name, ncores, bowtie2_index,
	cluster, jobs, keep_duplicates, max_javamem, extract_mito, reference_genome,
	clipl, clipr, keep_temp_files, skip_fastqc,
	bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
	bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path):
	
	"""
	proatac: a toolkit for prePROcessing ATAC-seq and scatac-seq data. \n
	Caleb Lareau; Buenrostro Lab. \n
	modes = ['bulk', 'check', 'single'] \n
	See http://proatac.readthedocs.io for more details.
	"""
	
	__version__ = get_distribution('proatac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting proatac pipeline v%s" % __version__)
	p = proatacProject(script_dir, mode, input, output, name, ncores, bowtie2_index,
		cluster, jobs, keep_duplicates, max_javamem, extract_mito, reference_genome,
		clipl, clipr, keep_temp_files, skip_fastqc,
		bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
		bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path)

	# -------------------------------
	# Atypical analysis modes
	# -------------------------------	
	
	if (mode == "check"):
		click.echo("Will process the following samples / files with bulk / single specified: \n\n")
		print("Sample", "Fastq1", "Fastq2")
		for x in range(len(p.samples)):
			print(p.samples[x], p.fastq1[x], p.fastq2[x])
		click.echo("\n\nIf this table doesn't look right, consider specifying a manually created sample input table.")
		sys.exit("Success!")

	# -------------
	# Adapter Trim
	# -------------
	trimfolder = outfolder + "/01_trimmed"
	if not os.path.exists(trimfolder + "_reads"):
		os.makedirs(trimfolder+ "_reads")
				
	click.echo(gettime() + "Trimming samples", logf)
	
	snakedict1 = {'allsamples' : parselfolder + "/allsamples.csv", 'outdir' : outfolder,
		'scriptdir' : script_dir}
		
	y1 = parselfolder + "/snake.trim.yaml"
	with open(y1, 'w') as yaml_file:
		yaml.dump(snakedict1, yaml_file, default_flow_style=False)
		
	snakecall1 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Trim --cores ' + p.max_cores + ' --config cfp="' + y1 + '"'
	os.system(snakecall1)
	click.echo(gettime() + "Sample trimming done.", logf)
	
	# -------------
	# Alignment
	# -------------
	if not os.path.exists(logfolder + "/bowtie2logs"):
		os.makedirs(logfolder+ "/bowtie2logs")
	if not os.path.exists(outfolder + "/02_aligned_reads"):
		os.makedirs(outfolder+ "/02_aligned_reads")
		
	click.echo(gettime() + "Aligning samples", logf)
	
	snakedict2 = {'bowtie2' : p.bowtie2_path, 'bowtie2index' : p.bowtie2_index,
		'outdir' : outfolder, 'samtools' : p.samtools_path}
	
	y2 = parselfolder + "/snake.align.yaml"
	with open(y2, 'w') as yaml_file:
		yaml.dump(snakedict2, yaml_file, default_flow_style=False)
	
	snakecall2 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.Align --cores ' + p.max_cores + ' --config cfp="' + y2 + '"'
	os.system(snakecall2)
	
	# -------------
	# bam process
	# -------------
	if not os.path.exists(outfolder + "/03_processed_reads"):
		os.makedirs(outfolder+ "/03_processed_reads")
	if not os.path.exists(outfolder + "/03_processed_reads/individual"):
		os.makedirs(outfolder+ "/03_processed_reads/individual")
	if not os.path.exists(logfolder + "/rmduplogs"):
		os.makedirs(logfolder+ "/rmduplogs")
		
	click.echo(gettime() + "Cleaning up .bam files", logf)
	
	# Determine chromosomes to keep / filter
	chrs = os.popen(p.samtools_path + " idxstats " +  outfolder + "/02_aligned_reads/* | cut -f1").read().strip().split("\n")
	rmchrlist = ["*", "chrY", "MT", "chrM"]
	keepchrs = [x for x in chrs if x not in rmchrlist and len(x) < int(p.chr_name_length)]
	
	snakedict3 = {'keepchrs' : keepchrs, 'read_quality' : p.read_quality, 'java' : p.java_path,
		'outdir' : outfolder, 'samtools' : p.samtools_path, 'project_name' : p.project_name, 'scriptdir' : script_dir}
		
	y3 = parselfolder + "/snake.bamprocess.yaml"
	with open(y3, 'w') as yaml_file:
		yaml.dump(snakedict3, yaml_file, default_flow_style=False)
	
	snakecall3 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.BamProcess --cores ' + p.max_cores + ' --config cfp="' + y3 + '"'
	os.system(snakecall3)
	
		
	# ---------------------
	# Peaks
	# ---------------------
	if not os.path.exists(outfolder + "/04_qc"):
		os.makedirs(outfolder+ "/04_qc")
		
	snakedict4 = {
		'bedtools' : p.bedtools_path, 'blacklistFile' : p.blacklistFile, 'macs2' : p.macs2_path,
		'macs2_genome_size' : p.macs2_genome_size, 'n_peaks' : p.n_peaks, 'outdir' : outfolder,
		'peak_width': p.peak_width, 'project_name' : p.project_name, 'R' : p.R_path, 'script_dir' : script_dir
	}	
		
	y4 = parselfolder + "/snake.callpeaksone.yaml"
	with open(y4, 'w') as yaml_file:
		yaml.dump(snakedict4, yaml_file, default_flow_style=False)
	
	snakecall4 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.CallPeaksOne --cores ' + p.max_cores + ' --config cfp="' + y4 + '"'
	os.system(snakecall4)
	
	# ---------------------
	# Individual QC
	# ---------------------	
	if not os.path.exists(outfolder + "/04_qc/individualQC"):
		os.makedirs(outfolder+ "/04_qc/individualQC")
		
	snakedict5 = {
		'bedtools' : p.bedtools_path, 'outdir' : outfolder,
		'project_name' : p.project_name, 'R' : p.R_path, 'samtools' : p.samtools_path, 
		'script_dir' : script_dir
	}	
		
	y5 = parselfolder + "/snake.qcstats.yaml"
	with open(y5, 'w') as yaml_file:
		yaml.dump(snakedict5, yaml_file, default_flow_style=False)
	
	snakecall5 = 'snakemake --snakefile ' + script_dir + '/bin/snake/Snakefile.QCstats --cores ' + p.max_cores + ' --config cfp="' + y5 + '"'
	os.system(snakecall5)
		
	# Suspend logging
	logf.close()
	
	
