#!/usr/bin/python3

#import the necessary modules
import NGSutilities
import os
import argparse
import shlex
import subprocess

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='runFastQC.py', description='Handle input arguments for calling FastQC')
	parser.add_argument('-i', '--input', help='Input file name or directory', required=True)
	parser.add_argument('-n', '--nthreads', help='Number of cores to be allocated. Defaults to 1 if not set')
	parser.add_argument('-g', '--genome', help='Location of the Star genome file. See the STAR manual for details', required=True)
	args = parser.parse_args()
	
	#get list of Fastq files
	fq = NGSutilities.findFiles(args.input, pattern='*.fastq*')
	
	cmd = 'star --genomeDir ' + args.genome + ' --runThreadN ' + args.nthreads + ' --outFileNamePrefix /home/data/pbt/RNASeq/BamFiles/SampledReads/subsetNMB261/ ' +	'--readFilesIn /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R1_001.fastq /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R2_001.fastq'

	args = shlex.split(cmd)
	done = subprocess.check_call(args)


#genome = '/home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/STARgenomes/hg19'
#nthreads = str(1)
#fq = NGSutilities.findFiles('/home/data/pbt/RNASeq/rawData/SampledReads', pattern='*.fastq*')
#cmd = 'star --genomeDir ' + genome + ' --runThreadN ' + nthreads + ' --outFileNamePrefix /home/data/pbt/RNASeq/BamFiles/SampledReads/subsetNMB261/ ' +	'--readFilesIn /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R1_001.fastq /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R2_001.fastq'

#commandline = shlex.split(cmd)
#done = subprocess.check_call(commandline)