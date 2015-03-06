# prioritize -b NBINS -t THRESHOLD GWAS_DATA -o DESTINATION

import argparse
import sys
import requests
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

# Check if the array is well-formed
# Need to check if coordinate is out-of-bounds too
def validate(arr):
	if len(arr)!=3:
		return False
	else:
		return True

def myProgressBar(maxval):
	widgets=[Percentage(), ' ', Bar()," ", AdaptiveETA(), ' ', Timer()]
	return ProgressBar(widgets=widgets, maxval=maxval)

### CONSTANTS
CANYON_API="http://localhost:3000/genome"

### Parsing Command Line Options
parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
	                          help="Path to GWAS Data")
parser.add_argument("-o", metavar="DESTINATION_PATH",
	                        default="result.data",
	                        help="Path to output file, default to result.data")
parser.add_argument("-b", metavar="NBINS", 
	                        type=int, help="Number of bins")
parser.add_argument("-t", metavar="THRESHOLD", type=float, default=0.5,
	                        help="Threshold, range in (0,1)")
args = parser.parse_args()
if args.t<=0 or args.t>=1:
	print("Prioritize: Invalid Threshold. Range should be (0,1)")
	sys.exit(1)

### Reading GWAS Data
print("Reading GWAS Data...")
with open(args.gwas,'r') as data_file:
	nlines = sum(1 for _ in data_file)
gwas_data = []
pbar = myProgressBar(nlines)
with open(args.gwas,'r') as data_file:
	for line in pbar(data_file):
		base = line.split('\t')
		if not validate(base):
			print("Prioritize: Invalid GWAS Data: {0}".format(line))
			sys.exit(1)
		else:
			gwas_data.append(Genome(int(base[0]), int(base[1]), Decimal(base[2])))

### Downloading Canyon Data
gwas_data.sort(key=lambda g: (g.chrom, g.loc))
canyon_data = []
print("Downloading Canyon Data...")
pbar = myProgressBar(23)
pbar.update(0)
for chrom, group in groupby(gwas_data, key=lambda g: g.chrom):
	locs = map(lambda g: (g.loc), group)
	res = requests.post(CANYON_API, timeout=120,
		                  data={'chrom':chrom, 'locs[]': locs})
	pbar.update(chrom)
	canyon_data.extend(map(Decimal, res.json()))
pbar.finish()

gwas_value = map(lambda g: g.value, gwas_data)


### Computation Starts Here with gwas_value and canyon_data
print("Writing Result to File...")
with open(args.o,'w') as output_file:
	for data1, data2 in zip(gwas_value,canyon_data):
		output_file.write("{0}\t{1}\n".format(data1,data2))



