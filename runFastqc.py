#!/usr/bin/python3
"""
Wrapper to the FastQC program available from http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Run the FastQC program on a list of Fastq files. This requires that FastQC is installed on the current system

"""

#import the necessary modules
import NGSutilities
import os
import argparse
import shlex
import subprocess

#function definition
def runFastQC(fqFiles, reportdir, nthread=1):
	"""Run the FastQC program on a list of Fastq files
	This requires that FastQC is installed on the current system
	
	Arguments
	---------
	fqFiles - List of Fastq files, either zipped or uncompressed
	
	reportdir - Directory path where FastQC report files are saved
	
	nthread - Number of threads FastQC can use to process files in parallel. Default=1
	"""
	
	#construct the fastqc command
	cmd = 'fastqc -t ' +  nthread + ' -o ' + reportdir + ' ' + ' '.join(fqFiles)
	
	args = shlex.split(cmd)
	qc = subprocess.check_call(args)
	
	print('QC reporting complete')

#set up the output directory
output = '/home/data/pbt/RNASeq/QC/FastQCreport'
if not os.path.exists(output):
	os.makedirs(output)

if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='runFastQC.py', description='Handle input arguments for calling FastQC')
	parser.add_argument('-i', '--input', help='Input file name or directory', required=True)
	parser.add_argument('-n', '--nthreads', help='Number of cores to be allocated. Defaults to 1 if not set')
	args = parser.parse_args()

	#get list of Fastq files, run FastQC on the fastq files
	fq = NGSutilities.findFiles(args.input, pattern='*.fastq*')
	runFastQC(fq, nthread=args.nthreads, reportdir=output)

	#remove the zip files that fastqc creates by default as we never use these
	print('Remove unused Zip files')
	zipfiles = NGSutilities.findFiles(output, pattern='*fastqc.zip')
	rmfiles= [os.remove(f) for f in zipfiles]
