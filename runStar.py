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
	
	#build output directory path. If it does not exist create it
	batchname = os.path.basename(args.input)
	output = os.path.join('/home', 'data', 'pbt', 'RNASeq', 'BamFiles', batchname)
	if not os.path.exists(output):
		os.makedirs(output)
	
	output = output + '/'
	
	#organise the list of files to be compatible with STAR
	#split this into 2 lists, left reads and right reads using conditional list comprehension
	left = [f for f in fq if 'R1' in f]
	right = [f for f in fq if 'R2' in f]
	
	#sort both lists. The order is essentially alphabetical. Important factor is that
	#the paired read files are in the same position in their respective lists
	left.sort()
	right.sort()
	
	#put them back together again as a comma separated string
	sortedFiles = ','.join(left + right)
	
	#construct the command
	alignment = 'star --genomeDir ' + args.genome + ' --runThreadN ' + args.nthreads + ' --outFileNamePrefix ' + output +	' --readFilesIn /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R1_001.fastq /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R2_001.fastq' + ' --outStd SAM'
	
	#run Star, send standard out to a pipe to go to samtools below
	align = shlex.split(alignment)
	sam = subprocess.Popen(align, stdout=subprocess.PIPE)
	
	bam = 'samtools view -bS -o ' + output + 'Sample_NMB261.bam ' + '-'
	bam = shlex.split(bam)
	unsortedBam = subprocess.Popen(bam, stdin=sam.stdout)
	sam.stdout.close()
	ret = unsortedBam.communicate()[0]


#genome = '/home/data/genomes/Homo_sapiens/UCSC/hg19/Sequence/STARgenomes/hg19'
#nthreads = str(1)
#fq = NGSutilities.findFiles('/home/data/pbt/RNASeq/rawData/SampledReads', pattern='*.fastq*')


#cmd = 'star --genomeDir ' + genome + ' --runThreadN ' + nthreads + ' --outFileNamePrefix /home/data/pbt/RNASeq/BamFiles/SampledReads/subsetNMB261/ ' +	'--readFilesIn /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R1_001.fastq /home/data/pbt/RNASeq/rawData/SampledReads/Sample_NMB261/subsetNMB261_TAGCTT_L004_R2_001.fastq'

#commandline = shlex.split(cmd)
#done = subprocess.check_call(commandline)