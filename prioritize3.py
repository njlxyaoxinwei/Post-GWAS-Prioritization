# prioritize [-b NBINS] [-t THRESHOLD] [-o DESTINATION] [-a ANNOTATION] GWAS_DATA 

import argparse
import sys
import requests
import numpy as np
from itertools import groupby
from decimal import *
from progressbar import *


class Genome:
	def __init__(self, chrom, loc, value):
		self.chrom = chrom
		self.loc = loc
		self.value = value
	def __str__(self):
		return "{0}:{1}:{2}".format(self.chrom, self.loc, self.value)

def line_count(path):
	with open(path,'r') as f:
		nlines = sum(1 for _ in f)
	return nlines

def die(string):
	print("Prioritize: " + string)
	exit(1)

# checks if coordinate is out-of-bounds
def validate(chrom, loc):
	if chrom in range(1,25):
		return True
	else:
		return False

def myProgressBar(maxval):
	widgets=[Percentage(), ' ', Bar()," ", AdaptiveETA(), ' ', Timer()]
	return ProgressBar(widgets=widgets, maxval=maxval)

def get_canyon(gwas_data):
	canyon_data = []
	pbar = myProgressBar(24)
	pbar.update(0)
	for chrom, group in groupby(gwas_data, key=lambda g: g.chrom):
		locs = map(lambda g: (g.loc), group)
		res = requests.post(CANYON_API, timeout=120,
		                  	data={'chrom':chrom, 'locs[]': locs})
		pbar.update(chrom)
		canyon_data.extend(map(Decimal, res.json()))
	pbar.finish()
	return canyon_data

def get_parser():
	parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
	parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
		                          help="Path to GWAS Data")
	parser.add_argument("-o", metavar="DESTINATION_PATH",
		                        default="result.data",
		                        help="Path to output file, default to result.data")
	parser.add_argument("-b", metavar="NBINS", 
		                        type=int, help="Number of bins")
	parser.add_argument("-t", metavar="THRESHOLD", type=float, default=0.1,
		                        help="Threshold, range in (0,1)")
	parser.add_argument("-a", metavar="ANNOTATION_PATH",
		                        help="Path to functional annotation file," 
		                             " default to data from GenoCanyon")
	return parser	

def parse_line(line):
	arr = line.split('\t')
	if len(arr)==3:
		chrom, loc, value = int(arr[0]), int(arr[1]), Decimal(arr[2])
		if validate(chrom, loc):
			return chrom, loc, value
	die("Invalid data: {0}".format(line))

def get_data(path):
	nlines = line_count(path)
	pbar = myProgressBar(nlines)
	data = []
	with open(path,'r') as f:
		for line in pbar(f):
			data.append(Genome(*parse_line(line)))
	data.sort(key=lambda g: (g.chrom, g.loc))
	return data



################################################################################
####                              Main Program                              ####
################################################################################

######  CONSTANTS 
#===============================================================================
CANYON_API="http://localhost:3000/genome"
#===============================================================================

######  INPUT HANDLING
#===============================================================================
### Parsing Command Line Options
args = get_parser().parse_args()

### Input Validation
if args.t<=0 or args.t>=1:
	die("Invalid THRESHOLD {0}. Range should be (0,1)".format(args.t))

if args.b!=None and args.b<=0:
	die("Invalid NBINS {0}. Should be a positive integer".format(args.b))

### Reading GWAS Data
print("Reading GWAS Data...")
gwas_data = get_data(args.gwas)

### Getting Canyon Data
if args.a==None:
	print("Downloading Canyon Data...")
	canyon_scores = get_canyon(gwas_data)
else:
	print("Reading Annotation Data...")
	canyon_scores = list(map(lambda g: g.value, get_data(args.a)))

### Prepare for computation
pvalues = list(map(lambda g: g.value, gwas_data))
count1, count2 = len(pvalues), len(canyon_scores)
if count1!=count2:
	die("Pvalues and Annotation scores do not match: "
		  "{0} pvalues and {1} scores".format(count1, count2))
pvalue_v = np.array(pvalues, dtype=np.float64)
canyon_v = np.array(canyon_scores, dtype=np.float64)

###### COMPUTATION
#===============================================================================
print(len(pvalue_v))
print(len(canyon_v))
exit(0)



### Computation Starts Here with gwas_value and canyon_data
print("Writing Result to File...")
with open(args.o,'w') as output_file:
	for data1, data2 in zip(gwas_value,canyon_data):
		output_file.write("{0}\t{1}\n".format(data1,data2))



