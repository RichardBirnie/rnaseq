#!/usr/bin/python3

#import the necessary modules
import NGSutilities
import os
import argparse
import shlex
import subprocess

#get list of Fastq files
fq = NGSutilities.findFiles(args.input, pattern='*.fastq*')
