#!/usr/bin/python3

""" A collection of functions to support the RNAseq processing pipeline"""

#import some useful modules
import fnmatch
import os

#define necessary functions
def findFiles(path, pattern='*.fastq.gz'):
	""" Recursively search all directories below 'path' that match 'pattern'.
		return the results as a list """
	
	#create a list of files to store output based on a search pattern
	Files = []
	for root, dirnames, filenames in os.walk(path):
		for filename in fnmatch.filter(filenames, pattern):
			Files.append(os.path.join(root, filename))
	
	return Files
