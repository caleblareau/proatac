import pysam
import click
import os
from pkg_resources import get_distribution

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['ATAC','RNA', 'PE', 'SE', 'atac', 'rna', 'pe', 'se']))

@click.option('--input', '-i', default = ".", required=True, help='.bam file to be converted to .bigwig file')


def main(mode, input):
	
	"""
	pybam2bw: Convert bam files into bigwigs with SPMR normalization \n
	Caleb Lareau
	"""
	
	__version__ = get_distribution('proatac').version
	input_filename = input
	click.echo("Starting bam to bigwig conversion")

	# Create file of contig lengths
	contigfile = input_filename.replace(".bam", ".chr.lengths.tsv")
	
	bam_input = pysam.AlignmentFile(input_filename, 'rb')
	header = str(bam_input.header).split("\n")
	contigfile_handle = open(contigfile, "w")

	for i in header:
		i_list = i.split("\t")
		if(i_list[0] == '@SQ'):
			out1 = i_list[1].replace("SN:", "")
			out2 = i_list[2].replace("LN:", "")
			contigfile_handle.write(out1 + "\t" + out2 + "\n")

	contigfile_handle.close()
	
	# Create bdg file
	click.echo("Running macs2 to assemble .bdg file")
	if(mode in ["ATAC", "atac", "PE", "pe"]):
		format = " -f BAMPE"
	else:
		format = " -f BAM"
	macs2call = "macs2 callpeak -t " + input_filename + " -B --SPMR -n " + input_filename + format
	print(macs2call)
	os.system(macs2call)	

	# Threshold and slop
	click.echo("Adjusting bdg")
	treat_pileup = input_filename + "_treat_pileup.bdg"
	bedtoolscmd = "bedtools slop -i "+treat_pileup+" -g "+contigfile+" -b 0 | bedClip stdin "+contigfile+" "+treat_pileup+".clip"
	os.system(bedtoolscmd)	
	
	# Sort
	sortcmd = "LC_COLLATE=C sort -k1,1 -k2,2n "+treat_pileup+".clip > "+treat_pileup+".sort.clip"
	os.system(sortcmd)	
	
	# Final bw
	finalbw = input_filename.replace(".bam", ".bw")
	bg2bwcmd = "bedGraphToBigWig " + treat_pileup+".sort.clip "+contigfile+ " " + finalbw
	os.system(bg2bwcmd)	
	
	# Remove temp
	#os.remove(treat_pileup)
	#os.remove(treat_pileup+".sort.clip")
	#os.remove(treat_pileup+".clip")
	#os.remove(treat_pileup)

	
	
	
	
	
	
	