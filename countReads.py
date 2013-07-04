#!/usr/bin/python3
"""
Assign reads to genes or exons using the HTSeq toolkit http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html 

This program takes the following basic steps:

Sort bam file into queryname order using samtools sort. This puts paired reads on adjacent lines. The output from this is piped directly to samtools view

samtools view converts the sorted bam to sam and pipes it directly to one of either htseq-count or dexseq_count which assigns reads to genes or exons respectively

Each bamfile is handled in a separate process so multiple files can be counted in parallel by setting the -n argument. Type python3 countReads.py -h for commandline usage instructions

"""


import NGSutilities
import os
import argparse
import shlex
import subprocess
from multiprocessing import Pool

#function definitions
def compileCommands(files):
	"""Takes a list of bamfiles and constructs a samtools command to sort each file in read name order ready for htseq-count and returns a list of these commands. See help(countReadsGenes.compileCommands) for usage. This is the starting point for branching multiple processes
	
	Arguments
	---------
	files - a list of bam files to be processed with full absolute file paths
	
	Return
	------
	A list of samtools sort command strings ready to be passed to runCommand()
	
	"""
	commands = []
	for f in files:
		namesort = 'samtools sort -no ' + f + ' ' + os.devnull
		commands.append(namesort)
	
	return commands

def countGenes(c):
	""" Takes a samtools sort command, executes it, and pipes the output to samtools view.
	
	samtools view converts the stream from bam to sam and pipes it through to htseq-count. 
	
	htseq-count assigns each read to a gene based on the contents of a gtf file describing the coordinates of genes and other features in the genome with the htseq-count mode option set to intersection-nonempty
	
	Arguments
	---------
	c - a single samtools sort command string
	
	Usage
	-----
	This function is more typically called as:
	p = Pool(processes=int(args.ncores))
	outputs = p.map(runCommand, commands)
	
	This creates a pool of worker processes and sends one command to each process.
	"""
	
	#sort the file by query name, i.e. read name.
	#This is required by htseq-count to detect paired reads
	#pipe output to samtools view to get it into SAM format
	namesort = shlex.split(c)
	infile = namesort[3]
	nsorted = subprocess.Popen(namesort, stdout=subprocess.PIPE)
	
	#catch stdout from samtools sort, convert to sam and pipe directly to htseq-count
	sam = 'samtools view -h - '
	sam = shlex.split(sam)
	sam =subprocess.Popen(sam, stdin=nsorted.stdout, stdout=subprocess.PIPE)
	nsorted.stdout.close()
	
	#generate output path. If the directory does not exist create it
	output = infile.replace('BamFiles', 'ReadCounts/Genes')
	output = output.replace('sortedRG.bam', 'geneCounts.txt')
	output = os.path.join(os.path.dirname(os.path.dirname(output)), os.path.basename(output))
	if not os.path.exists(os.path.dirname(output)):
		os.makedirs(os.path.dirname(output), exist_ok=True)
	
	#get sam from pipe above and count reads
	countreads = 'htseq-count -s ' + args.stranded + ' -t exon -m intersection-nonempty - ' + args.gtf_file
	countreads = shlex.split(countreads)
	cts = subprocess.Popen(countreads, stdin=sam.stdout, stdout=open(output, 'w'))
	sam.stdout.close()
	ret = cts.communicate()

def countExons(c):
	""" Takes a samtools sort command, executes it, and pipes the output to samtools view.
	
	samtools view converts the stream from bam to sam and pipes it through to dexseq_count. 
	
	dexseq_count assigns each read to an exon/counting bin. This requires a modified gtf file that has been processed with the tool dexseq_prepare_annotation.py from the DEXSeq Bioconductor package to split overlapping exons into counting bins. The initial gtf prior to modification should describe the coordinates of genes, exons and other annotations in the genome such as those available from the GENCODE project http://www.gencodegenes.org/
	
	Arguments
	---------
	c - a single samtools sort command string
	
	Usage
	-----
	This function is more typically called as:
	p = Pool(processes=int(args.ncores))
	outputs = p.map(runCommand, commands)
	
	This creates a pool of worker processes and sends one command to each process.
	"""
	#sort the file by query name, i.e. read name.
	#This is required by htseq-count to detect paired reads
	#pipe output to samtools view to get it into SAM format
	namesort = shlex.split(c)
	infile = namesort[3]
	nsorted = subprocess.Popen(namesort, stdout=subprocess.PIPE)
	
	#catch stdout from samtools sort, convert to sam and pipe directly to htseq-count
	sam = 'samtools view -h - '
	sam = shlex.split(sam)
	sam =subprocess.Popen(sam, stdin=nsorted.stdout, stdout=subprocess.PIPE)
	nsorted.stdout.close()
	
	#generate output path. If the directory does not exist create it
	output = infile.replace('BamFiles', 'ReadCounts/Exons')
	output = output.replace('sortedRG.bam', 'exonCounts.txt')
	output = os.path.join(os.path.dirname(os.path.dirname(output)), os.path.basename(output))
	if not os.path.exists(os.path.dirname(output)):
		os.makedirs(os.path.dirname(output))
	
	#get sam from pipe above and count reads
	cts = 'dexseq_count -s ' + args.stranded + ' -p yes ' + args.gtf_flat + ' - ' + output
	cts = shlex.split(cts)
	#call the command. Not the redirection of stdout and stderr to /dev/null. 
	#This command generates lots of output (especially in parallel) which blocks the stdout and stderr buffer eventually
	#This has the side effect of supressing all messages
	cts = subprocess.Popen(cts, stdin=sam.stdout, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
	#cts = subprocess.Popen(cts, stdin=sam.stdout)
	sam.stdout.close()
	ret = cts.communicate()

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='countReads.py', description='Use the HTseq tools to assign reads per gene and per exon for a set of bam files. http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html. This requires that the HTSeq toolkit has been installed on the current system and can be accessed by typing htseq-count at the commandline. In addition the script dexseq_count from the Bioconductor package DEXSeq must also be installed and accessible by typing dexseq_count at the commandline')
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated. Defaults to 1 if not set', default=1)
	
	parser.add_argument('-g', '--gtf_file', help='gtf file describing the coordinates of genes and other features in the genome', required=True)
	
	parser.add_argument('-f', '--gtf_flat', help="'Flattened' version of the gtf file used for the -g argument. This can be prepared with the dexseq_prepare_annotation.py script shipped as part of the DEXSeq Bioconductor package. Briefly, that script breaks up overlapping exons into discrete non-overlapping counting bins", required=True)
	
	parser.add_argument('-s', '--stranded', help="Allows setting of the htseq-count -s flag to specify whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: no). 'reverse' means 'yes' with reversed strand interpretation", default='no')
	
	args = parser.parse_args()
	
	#sanity check input file name
	#remove trailing / if there is one
	args.input = args.input.rstrip('/')
	
	#get batch name from the input directory
	batchname = os.path.basename(args.input)
	
	#get list of input files
	bf = NGSutilities.findFiles(args.input, pattern='*RG.bam')
	commands = compileCommands(bf)
	
	#Use the map function of the Pool object, not the default function
	#this distributes the jobs of map() across the processes available
	#to the Pool.
	#As with the standard implementation of map(), it takes a function and 
	#a list as arguments, and executes the function with every member of the
	#list as the sole arguement.
	print('Starting Counting Reads per Gene')
	#make a Pool object from the multiprocessing package.
	#processes argument specifies how many cores to use
	p1 = Pool(processes=int(args.ncores))
	outputs = p1.map(countGenes, commands)
	print('Finished Counting Reads per Genes')
	
	print('Starting Counting Reads per Exon')
	#make a Pool object from the multiprocessing package.
	#processes argument specifies how many cores to use
	p2 = Pool(processes=int(args.ncores))
	outputs = p2.map(countExons, commands)
	print('Finished Counting Reads per Exon')
