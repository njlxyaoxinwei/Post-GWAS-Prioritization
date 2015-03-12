#prioritize [-b NBINS] [-t THRESHOLD] [-o DESTINATION] [-a ANNOTATION] GWAS_DATA 
from __future__ import division
import argparse
import sys
import time
import requests
import numpy as np
from scipy import special
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

def make_request(data, times):
	try:
		res = requests.post(CANYON_API, timeout=30*times, data=data)
		res.raise_for_status()
		return res.json()
	except Exception as inst:
		print(inst)
		if times==MAX_RETRY:
			print("{0} attempts failed, ".format(MAX_RETRY) + 
				    "Please check Internet Connection")
			die("Network Error")
		else:
			print("Attempting again in 5 seconds...")
			time.sleep(5)
			return make_request(data, times+1)

def get_canyon(gwas_data):
	canyon_data = []
	pbar = myProgressBar(24)
	pbar.update(0)
	for chrom, group in groupby(gwas_data, key=lambda g: g.chrom):
		locs = list(map(lambda g: (g.loc), group))
		# CANNOT use map here otherwise only the first attempt can go through,
		# since the iterator will have finished the whole pass after that.
		res = make_request({'chrom': chrom, 'locs[]': locs},1)
		pbar.update(chrom)
		canyon_data.extend(map(Decimal, res))
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

def cross_validate_nbins(vector):
	def risk_func(hist, w, n):
		hist = np.array(hist)
		return np.sum(hist**2) / (n**2*w) - (2/(n*(n-1)*w)) * np.sum(hist*(hist-1))

	max_nbins = MAX_NBINS
	n = len(vector)
	a,b = 0,1 # Lower and Upper bound for data
	span = b-a
	risk = []
	pbar = myProgressBar(max_nbins)
	nbins_v = np.arange(max_nbins+1)[1:]

	for candidate in pbar(nbins_v):
		bins = np.linspace(a,b,candidate+1)
		hist = my_histogram(vector, nbins=candidate, density=False)
		risk.append(risk_func(hist, span/candidate, n))

	risk = np.array(risk)
	# Taking the largest nbin with the minimal risk
	result = nbins_v[risk==risk.min()][-1]
	return result

def get_bin(value_v, nbins, a=0, b=1):
	index_v = np.empty_like(value_v, dtype=int)
	np.ceil((value_v-a)/(b-a)*nbins, index_v)
	index_v = index_v - 1
	return index_v

def my_histogram(data_v, nbins, density=False, a=0, b=1):
	index_v = get_bin(data_v, nbins, a, b)
	hist = np.zeros(nbins, dtype=int)
	for index in index_v:
		hist[index]+=1
	if density:
		result = np.empty(nbins, dtype=FLOAT_TYPE)
		np.divide(hist*nbins, len(data_v)*(b-a), result)
		return result
	else:
		return hist
			
def conditional_exp(data_v, th0, th1, density_v, nbins):
	part1 = np.power(data_v, th1-1)/special.beta(th1,1)
	part2 = density_v[get_bin(data_v, nbins)]
	return th0*part1/(th0*part1+(1-th0)*part2)

def run_EM(pvalue_func_v, density_non_v, nbins):
	theta = np.array([0.01, 0.5], dtype=FLOAT_TYPE)
	theta_trace = []
	pbar = myProgressBar(MAX_ITER)
	for _ in pbar(range(MAX_ITER)):
		t1 = conditional_exp(pvalue_func_v, theta[0], theta[1], density_non_v, nbins)
		new_theta = np.empty(2, dtype=FLOAT_TYPE)
		new_theta[0] = np.mean(t1)
		new_theta[1] = new_theta[0]/np.mean(t1*(-np.log(pvalue_func_v)))
		theta = new_theta
		theta_trace.append(new_theta)
	return {'theta': theta, 'trace': theta_trace}




################################################################################
####                              Main Program                              ####
################################################################################
######  CONSTANTS 
#===============================================================================
CANYON_API="http://localhost:3000/genome"
MAX_RETRY=3
MAX_NBINS=100
FLOAT_TYPE = np.float64
MAX_ITER=5000
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
pvalue_v = np.array(pvalues, dtype=FLOAT_TYPE) #Pvalue vector
canyon_v = np.array(canyon_scores, dtype=FLOAT_TYPE) #Canyon score vector
print("Will perform Prioritization on {0} values".format(count1))
thd = args.t # Threshold for functional/non-functional divide
print("Will use {0} as threshold.".format(thd))
nbins = args.b
#===============================================================================
###### COMPUTATION
#===============================================================================
pvalue_func_v = pvalue_v[canyon_v > thd] # Functional Pvalues
pvalue_non_v = pvalue_v[canyon_v <= thd] # Non-functional Pvalues

### Cross-Validation
if nbins==None:
	print("Cross Validation...")
	nbins = cross_validate_nbins(pvalue_non_v)
	print("Will use {0} bins for optimal result".format(nbins))
else:
	print("Will use {0} bins as specified".format(nbins))
density_non_v = my_histogram(pvalue_non_v, nbins, density=True)

### EM
print("Running EM algorithm...")
result_EM = run_EM(pvalue_func_v, density_non_v, nbins)
[th0, th1] = result_EM['theta']
print(th0, th1)
### Posterior Calculation
# print("Calculating final result...")
# part1 = np.power(pvalue_v, (th1-1))/special.beta(th1, 1)
# part2 = density_non_v[get_bin(pvalue_v, nbins)]
# prior = th0*canyon_v
# posterior = prior*part1/(prior*part1 + (1-prior)*part2)


exit(0)

