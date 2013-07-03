#!/usr/bin/python3
"""
Wrapper to the STAR alignment program available from http://code.google.com/p/rna-star/

Run the STAR program on a list of Fastq files. This requires that STAR is installed on the current system and can be accessed via typing star at the commandline. This can be done by creating a symlink to the STAR install on linux systems something like 'sudo ln -s /path/to/star/installation /usr/local/bin/star'

This script can also be imported as a module by typing import runStar at the python interpreter provided this script is in your PYTHONPATH. This will give you access to the internal functions for interactive work

For commandline usage instructions type python3 runStar.py -h

Functions
---------
compileCommands - Compile a list of command strings to call STAR with appropriate settings on paired end reads for a single sample. See help(runStar.compileCommands) for usage.

runCommand - This function takes a single string describing a complete STAR command for 1 pair of sequencing read files corresponding to a single sample. The function will run the star command then pipe the output to samtools which converts sam to bam then sorts the bam in read name order.

"""

#import the necessary modules
import NGSutilities
import os
import argparse
import shlex
import subprocess
from multiprocessing import Pool

#function definitions
def compileCommands(files, genome, output):
	"""Compile a list of command strings to call STAR with appropriate settings on paired end reads for a single sample.
	
	Arguments
	---------
	files - a list of fastq files to be processed
	
	genome - string giving the directry path to the STAR genome file. See the STAR website for details
	
	output - Prefix to prepend at the start of the filename. Used to set the star --outFileNamePrefix argument
	"""
	#split the input list of files into 2 lists, left reads and right reads using conditional list comprehension
	left = [f for f in files if 'R1' in f]
	right = [f for f in files if 'R2' in f]
	
	#sort both lists. The order is essentially alphabetical. Important factor is that
	#the paired read files are in the same position in their respective lists
	left.sort()
	right.sort()
	
	commands = []
	for i in range(len(left)):
		#create output path
		n = os.path.basename(os.path.dirname(left[i]))
		savepath = os.path.join(output, n, n)
		
		#create the output directory
		if not os.path.exists(os.path.dirname(savepath)):
			os.makedirs(os.path.dirname(savepath))
		
		#construct the command
		alignment = 'star --genomeDir ' + genome + ' --runThreadN 1 --outFileNamePrefix ' + savepath + ' --readFilesIn ' + left[i] + ' ' + right[i] + ' --outStd SAM --outFilterMismatchNmax 2 --genomeLoad LoadAndKeep --outSAMstrandField intronMotif'
		commands.append(alignment)
		
	return commands

def runCommand(c):
	"""Runs a single shell command in a new subprocess
	
	This function takes a single string describing a complete STAR command for 1 pair of sequencing read files corresponding to a single sample. The function will run the star command then pipe the output to samtools which converts sam to bam then sorts the bam in read name order"""
	
	#run Star, send standard out to a pipe to go to samtools below
	align = shlex.split(c)
	samplename = align[9]
	sam = subprocess.Popen(align, stdout=subprocess.PIPE)
	
	#parse sample name to control output
	d = os.path.dirname(samplename)
	d = d.replace('rawData', 'BamFiles')
	output = os.path.join(d, os.path.basename(d))
	
	if not os.path.exists(os.path.dirname(output)):
		os.makedirs(os.path.dirname(output))
	
	#construct and call the samtools command to convert sam to bam
	bam = 'samtools view -bS -o ' + output + '.bam' + ' -'
	bam = shlex.split(bam)
	unsortedBam = subprocess.Popen(bam, stdin=sam.stdout)
	sam.stdout.close()
	ret = unsortedBam.communicate()[0]
	
	##THIS needs to be rearranged. Send unsorted Bam file directly to AddOrReplaceReadGroups and have this output in coordinate order
	##index the bam file in coordinate order
	##Alternatively don't do any sorting at this stage. Just add the read groups and validate the file format
	##sorting and indexing can be handled in runAlignmentQC to satisfy RNASeQC which appears to be the most pedantic
	
	##sort the bam file in coordinate order ready for RNASeQC
	##remove the unsorted file
	##sortedBamFile = output + '_sorted'
	##sortBam = 'samtools sort ' + output + ' ' + sortedBamFile
	##print('Sorting Bam file: ' + sortBam)
	##sortBam = shlex.split(sortBam)
	##sortedBam = subprocess.check_call(sortBam)
	##os.remove(output)
	
	#add read groups prior to validation
	#Nte that this sorts the bam file in coordinate order on the way out
	sampleID = os.path.basename(output)
	rg = os.path.basename(os.path.dirname(os.path.dirname(output)))
	addRG = 'picard-tools AddOrReplaceReadGroups INPUT=' + output + '.bam' + ' OUTPUT=' + output + '_sortedRG.bam' + ' SORT_ORDER=coordinate QUIET=TRUE RGID=' + rg + ' RGLB=Unknown RGPL=illumina RGPU=Unknown RGSM=' + sampleID
	print('Adding Read Groups: ' + output + '.bam')
	addRG = shlex.split(addRG)
	addedRG = subprocess.check_call(addRG)
	os.remove(output + '.bam')
	
	#validate bam file to check for consistency with the sam format specification
	validateBam = 'picard-tools ValidateSamFile INPUT=' + output + '_sortedRG.bam' + ' OUTPUT=/dev/stdout' + ' IGNORE=MISSING_TAG_NM QUIET=TRUE'
	print('Validating Bam file:' + output + '_sortedRG.bam')
	validateBam = shlex.split(validateBam)
	validatedBam = subprocess.Popen(validateBam, stdout=subprocess.PIPE)
	validationResult = validatedBam.communicate()[0]
	validationResult = str(validationResult, encoding='utf8')
	
	#check contents of validation test and report accordingly
	if validationResult.rstrip('\n') == 'No errors found':
		print(output + '_sortedRG.bam' + ': Validation check OK')
	else:
		print('VALIDATION ERROR: ' + output + '_sortedRG.bam' + '\n' + validationResult)


if __name__ == "__main__":
	#parse commandline arguments
	parser = argparse.ArgumentParser(prog='runStar.py', description='Use STAR to align a set of paired end fastq files. The resulting bam files are indexed using samtools index and have read groups added using picard AddOrReplaceReadGroups. The final bam file is sorted in genome coordinate order and validated against the sam format specification using picard ValidateSamFile.')
	
	parser.add_argument('-i', '--input', help='Input directory. Do not include the trailing / at the end of the filename as this will affect the internal operation of the script', required=True)
	parser.add_argument('-n', '--ncores', help='Number of cores to be allocated', required=True)
	parser.add_argument('-g', '--genome', help='Location of the Star genome file. See the STAR manual for details', required=True)
	args = parser.parse_args()
	
	#get list of Fastq files
	fq = NGSutilities.findFiles(args.input, pattern='*.fastq*')
	
	#build output directory path. If it does not exist create it
	batchname = os.path.basename(args.input)
	output = os.path.join('/home', 'data', 'pbt', 'RNASeq', 'BamFiles', batchname)
	if not os.path.exists(output):
		os.makedirs(output)
	
	#compile commands
	commands = compileCommands(fq, genome=args.genome, output=output)
	
	#make a Pool object from the multiprocessing package.
	#processes argument specifies how many cores to use
	p = Pool(processes=int(args.ncores))
	
	#Use the map function of the Pool object, not the default function
	#this distributes the jobs of map() across the processes available
	#to the Pool.
	#As with the standard implementation of map(), it takes a function and 
	#a list as arguments, and executes the function with every member of the
	#list as the sole arguement.
	print('Starting Alignments')
	outputs = p.map(runCommand, commands)
	print('Alignments Complete')