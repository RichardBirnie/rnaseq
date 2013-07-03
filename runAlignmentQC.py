#!/usr/bin/python3

"""
Wrapper to the RNASeQC program available from http://www.broadinstitute.org/cancer/cga/rna-seqc

Run the RNASeQC program on a list of Bam files. 

This script can also be imported as a module by typing import runAlignmentQC at the python interpreter provided this script is in your PYTHONPATH. This will give you access to the internal functions for interactive work

For commandline usage instructions type python3 runAlignmentQC.py -h

Functions
---------
compileCommands -  Takes a list of bamfiles and constructs a samtools command to count the number of reads in each file and returns a list of these commands. See help(runAlignmentQC.compileCommands) for usage. This is the starting point for branching multiple processes

runCommand - This function takes a single string describing a samtools command to count reads in the bam file. This function will then extract whichever is smaller, either 1million reads or all reads for the specified bam file. The output is piped directly to ReorderSam which creates a temporary bam file sorted in lexographic order (chr1, chr10, chr11, chr12 etc) as required by RNASeQC. The temporary file is indexed using samtools index. The temporary file is then passed to RNASeQC. Once RNASeQC completes all temporary files are deleted.

This function is typically called as:
p = Pool(processes=int(args.ncores))
outputs = p.map(runCommand, commands)

This creates a pool of worker processes and sends one command to each process.

"""

#import the necessary modules
import NGSutilities
import os
import argparse
import shlex
import subprocess
from multiprocessing import Pool

#function definitions
def compileCommands(files):
	"""Takes a list of bamfiles and constructs a samtools command to count the number of reads in each file and returns a list of these commands. See help(runAlignmentQC.compileCommands) for usage. This is the starting point for branching multiple processes
	
	Arguments
	---------
	files - a list of bam files to be processed with full absolute file paths
	
	Return
	------
	A list of samtools command strings ready to be passed to runCommand()
	
	"""
	commands = []
	for f in files:
		countreads = 'samtools view -h -c ' + f
		commands.append(countreads)
	
	return commands

def runCommand(c):
	"""This function takes a single string describing a samtools command to count reads in the bam file. This function will then extract whichever is smaller, either 1million reads or all reads for the specified bam file. The output is piped directly to ReorderSam which creates a temporary bam file sorted in lexographic order (chr1, chr10, chr11, chr12 etc) as required by RNASeQC. The temporary file is indexed using samtools index. The temporary file is then passed to RNASeQC. Once RNASeQC completes all temporary files are deleted.
	
	Arguments
	---------
	c - a single samtools command string of the form samtools view -h -c filename
	
	Usage
	-----
	This function is typically called as:
	p = Pool(processes=int(args.ncores))
	outputs = p.map(runCommand, commands)

	This creates a pool of worker processes and sends one command to each process.
"""
	countreads = shlex.split(c)
	infile = countreads[-1]
	nreads = subprocess.Popen(countreads, stdout=subprocess.PIPE)
	#retrieve the count from stdout
	nreads = nreads.communicate()[0]
	
	#if <1million present use all of them else
	#calculate the % corresponding to 1 million reads
	if int(nreads) > 1000000:
		samplesize = 1000000 / int(nreads)
	else:
		samplesize = -1
	
	#downsample the bam file sending the output to stdout which
	#will go to picard ReorderSam below
	downsample = 'samtools view -bh -s ' + str(samplesize) + ' ' + infile
	downsample = shlex.split(downsample)
	ds = subprocess.Popen(downsample, stdout=subprocess.PIPE)
	
	#order bam file in lexographical order, i.e. chr1, chr10, chr11 etc
	#input on stdin coming from downsampling above
	sample = os.path.basename(os.path.dirname(infile))
	output = os.path.dirname(infile)
	savepath = os.path.join(output, sample) + '_tmp.bam'
	orderBam = 'picard-tools ReorderSam INPUT=/dev/stdin' + ' OUTPUT=' + savepath + ' REFERENCE=' + args.reference
	#orderBam = 'picard-tools ReorderSam INPUT=/dev/stdin' + ' OUTPUT=' + savepath + ' REFERENCE=' + reference
	orderBam = shlex.split(orderBam)
	orderedBam = subprocess.Popen(orderBam, stdin=ds.stdout)
	ds.stdout.close()
	ret = orderedBam.communicate()[0]
	
	#index the bam file
	indexBam = 'samtools index ' + savepath
	print('Index Bam file: ' + indexBam)
	indexBam = shlex.split(indexBam)
	indexedBam = subprocess.check_call(indexBam)
	
	#run QC
	#create output path
	d = os.path.dirname(infile)
	d = d.replace('BamFiles', 'QC/BamQC')
	if not os.path.exists(os.path.dirname(d)):
		os.makedirs(os.path.dirname(d))
	
	#specify locations of the annotation files describing coordinates of genes in the genome
	gtfFile = '/home/data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/gencode.v7.annotation.gtf'
	#compile string for the -s argument. This tells RNASeQC where to look for the bam file
	#quote pattern at start and end of line is needed to enclose the main string in double quotes
	s = '"' + '|'.join([sample, savepath, batchname]) + '"'
	runQC = 'java -jar /opt/RNA-SeQC/RNA-SeQC_v1.1.7.jar'  + ' -s ' + s + ' -o ' + d + ' -n 1000 -t ' + gtfFile + ' -r /home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'
	runQC = shlex.split(runQC)
	qc = subprocess.check_call(runQC)
	
	#remove the temporary sorted file
	tmp = [savepath, savepath + '.bai']
	ret = [os.remove(t) for t in tmp]

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='runAlignmentQC.py', description='Use the Broad Institute tool RNASeQC (http://www.broadinstitute.org/cancer/cga/rna-seqc) to generate quality control metrics for a series of Bam files')
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated. Defaults to 1 if not set', default=1)
	parser.add_argument('-r', '--reference', help='Reference sequence in fasta format. Used by picard-tools ReorderSam to prepare the bam file prior to running RNASeQC. A sequence dictionary corresponding to the reference sequence should be presenting in the same directory. See Picard manual for details ', required=True)
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
	print('Starting QC')
	outputs = p.map(runCommand, commands)
	print('QC Complete')