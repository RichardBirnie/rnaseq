#!/usr/bin/python3

import NGSutilities
import os
import argparse
import shlex
import subprocess
from multiprocessing import Pool

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='countReadsGenes.py', description='Use the HTseq script htseq-count to assign reads to genes for a set of bam files. http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html. This requires that the HTSeq toolkit has been installed on the current system and can be accessed by typing htseq-count at the commandline')
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated. Defaults to 1 if not set', default=1)
	
	parser.add_argument('-gtf', '--gtf_file', help='gtf file describing the coordinates of genes and other features in the genome', required=True)
	
	args = parser.parse_args()
	
	#get batch name from the input directory
	batchname = os.path.basename(args.input)
	
	#get list of input files
	bf = NGSutilities.findFiles(args.input, pattern='*RG.bam')
	
	#sort the file by query name, i.e. read name.
	#This is required by htseq-count to detect paired reads
	#pipe output to samtools view to get it into SAM format
	namesort = 'samtools sort -no ' + bf[0] + ' ' + os.devnull
	namesort = shlex.split(namesort)
	nsorted = subprocess.Popen(namesort, stdout=subprocess.PIPE)

	#catch stdout from samtools sort, convert to sam and pipe directly to htseq-count
	sam = 'samtools view -h - '
	sam = shlex.split(sam)
	sam =subprocess.Popen(sam, stdin=nsorted.stdout, stdout=subprocess.PIPE)
	nsorted.stdout.close()

	#generate output path. If the directory does not exist create it
	output = bf[0].replace('BamFiles', 'ReadCounts/Genes')
	output = output.replace('sortedRG.bam', 'geneCounts.txt')
	output = os.path.join(os.path.dirname(os.path.dirname(output)), os.path.basename(output))
	if not os.path.exists(os.path.dirname(output)):
		os.makedirs(os.path.dirname(output))

	#get sam from pipe above and count reads
	countreads = 'htseq-count -s no -t exon -m intersection-nonempty - ' + args.gtf_file
	#countreads = 'htseq-count -s no -t exon -m intersection-nonempty - ' + '/home/data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
	countreads = shlex.split(countreads)
	cts = subprocess.Popen(countreads, stdin=sam.stdout, stdout=open(output, 'w'))
	sam.stdout.close()
	ret = cts.communicate()
	
	
	