#!/usr/bin/python3

#import necessary modules
import os
import argparse
import shlex
import subprocess
import logging
import datetime

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='runPipeline.py', formatter_class=argparse.RawDescriptionHelpFormatter, description='''Top level control script to run the RNASeq processing pipeline on a collection of Fastq files.
	
	This script takes two arguments an input directory and the number of processor cores you wish to allocate.
	Each Fastq file then goes through the following steps:
	
	runFastqc - assess the quality of the raw data using the FastQC program 
	
	runStar - align the raw data to the reference human genome using the STAR aligner. The output is a bam file
	sorted in genome coordinate order and a corresponding index filename 
	
	runAlignmentQC - Assess the quality of the aligned data using the RNASeQC program. For efficiency if your data
	contains more then 1 million reads this module will sample 1 million random reads using samtools and run the QC on that
	sample. This is usually representative of the whole file with the exception that coverage calculations will
	have to be scaled accordingly 
	
	countReads - assign reads to genes or exons using either htseq-count or dexseq_count. Each input file produces
	two output files, one file of read counts per gene and one file of read counts per exon. Each row in these
	files is labelled with an Ensembl ID corresponding to the gene or exon. All modules work in parallel if the -n
	flag is set >1. Output files are created in subdirectories alongside the raw data directory set with the -i
	flag.'''
	)
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated. Defaults to 1 if not set', default=1)
	
	args = parser.parse_args()
	
	#setup logging
	logfile = os.path.join(os.path.dirname(os.path.dirname(args.input)), 'PipelineLog_' + datetime.date.today().isoformat() + '.log')
	logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO)
	
	print('################# Start pipeline ################')
	logging.info(datetime.datetime.now().ctime() + ': Start Pipeline')
	#run fastqc
	print('Run Fastqc')
	fastqc = 'python3 runFastqc.py -i ' + args.input + ' -n ' + args.ncores
	logging.info(datetime.datetime.now().ctime() + ': ' + fastqc) #log the command string
	fastqc = shlex.split(fastqc)
	fastqc = subprocess.check_call(fastqc)
	
	#run star alignment
	print('Run Star alignment')
	star = 'python3 runStar.py -i ' + args.input + ' -g /home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/STARgenomes/hg19 -n ' + args.ncores
	logging.info(datetime.datetime.now().ctime() + ': ' + star)
	star = shlex.split(star)
	star = subprocess.check_call(star)
	
	#run QC on bam file
	print('Run post-alignment QC')
	bam = args.input
	bam = bam.replace('rawData', 'BamFiles')
	bamqc = 'python3 runAlignmentQC.py -i ' +  bam + ' -r /home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -n ' + args.ncores
	logging.info(datetime.datetime.now().ctime() + ': ' + bamqc)
	bamqc = shlex.split(bamqc)
	bamqc = subprocess.check_call(bamqc)
	
	#count reads per gene and per exon
	print('Count Reads')
	count = 'python3 countReads.py -i ' + bam + ' -g /home/data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.gtf -f /home/data/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/Gencode17/gencode.v17.annotation.flattened.for.dexseq.gtf -n ' + args.ncores
	logging.info(datetime.datetime.now().ctime() + ': ' + count)
	count = shlex.split(count)
	count = subprocess.check_call(count)
	
	print('Pipeline complete. Output files can be found at ' + os.path.dirname(os.path.dirname(args.input)))
	logging.info(' Pipeline complete. Output files can be found at ' + os.path.dirname(os.path.dirname(args.input)))
	