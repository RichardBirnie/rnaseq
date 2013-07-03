#!/usr/bin/python3

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

def runCommand(c):
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
		os.makedirs(os.path.dirname(output))
	
	#get sam from pipe above and count reads
	countreads = 'htseq-count -s ' + args.stranded + ' -t exon -m intersection-nonempty - ' + args.gtf_file
	countreads = shlex.split(countreads)
	cts = subprocess.Popen(countreads, stdin=sam.stdout, stdout=open(output, 'w'))
	sam.stdout.close()
	ret = cts.communicate()

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='countReadsGenes.py', description='Use the HTseq script htseq-count to assign reads to genes for a set of bam files. http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html. This requires that the HTSeq toolkit has been installed on the current system and can be accessed by typing htseq-count at the commandline')
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated. Defaults to 1 if not set', default=1)
	
	parser.add_argument('-gtf', '--gtf_file', help='gtf file describing the coordinates of genes and other features in the genome', required=True)
	
	parser.add_argument('-s', '--stranded', help="Allows setting of the htseq-count -s flag to specify whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: no). 'reverse' means 'yes' with reversed strand interpretation", default='no')
	
	args = parser.parse_args()
	
	#get batch name from the input directory
	batchname = os.path.basename(args.input)
	
	#get list of input files
	bf = NGSutilities.findFiles(args.input, pattern='*RG.bam')
	commands = compileCommands(bf)
	
	#make a Pool object from the multiprocessing package.
	#processes argument specifies how many cores to use
	p = Pool(processes=int(args.ncores))
	
	#Use the map function of the Pool object, not the default function
	#this distributes the jobs of map() across the processes available
	#to the Pool.
	#As with the standard implementation of map(), it takes a function and 
	#a list as arguments, and executes the function with every member of the
	#list as the sole arguement.
	print('Starting Counting Reads')
	outputs = p.map(runCommand, commands)
	print('Finished Counting Reads')
