# prioritize -b NBINS -t THRESHOLD GWAS_DATA -o DESTINATION

import argparse
import sys
import requests
from itertools import groupby
from decimal import *

class Genome:
	def __init__(self, chrom, loc, value):
		self.chrom = chrom
		self.loc = loc
		self.value = value
	def __str__(self):
		return "{0}:{1}:{2}".format(self.chrom, self.loc, self.value)

def validate(arr):
	if len(arr)!=3:
		return False
	else:
		return True

### CONSTANTS
CANYON_API="http://localhost:3000/genome"

### Parsing Command Line Options
parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
	                          type=argparse.FileType('r'), 
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
data_file = args.gwas
gwas_data = []
for line in data_file:
	base = line.split('\t')
	if not validate(base):
		print("Prioritize: Invalid GWAS Data: {0}".format(line))
		sys.exit(1)
	else:
		gwas_data.append(Genome(int(base[0]), int(base[1]), Decimal(base[2])))
data_file.close()

### Downloading Canyon Data
gwas_data.sort(key=lambda g: (g.chrom, g.loc))
canyon_data = []
for chrom, group in groupby(gwas_data, key=lambda g: g.chrom):
	locs = map(lambda g: (g.loc), group)
	res = requests.post(CANYON_API, timeout=100,
		                  data={'chrom':chrom, 'locs[]': locs})
	canyon_data.extend(map(Decimal, res.json()))
gwas_data = map(lambda g: g.value, gwas_data)


### 
with open(args.o,'w') as output_file:
	for data1, data2 in zip(gwas_data,canyon_data):
		output_file.write("{0}\t{1}\n".format(data1,data2))



