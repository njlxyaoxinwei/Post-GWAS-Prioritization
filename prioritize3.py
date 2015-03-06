# prioritize -b NBINS -t THRESHOLD GWAS_DATA -o DESTINATION

import argparse

parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
	                          type=argparse.FileType('r'), 
	                          help="Path to GWAS Data")
parser.add_argument("-o", metavar="DESTINATION_PATH",
	                        type=argparse.FileType('w'),
	                        default=open("result.data",'w'),
	                        help="Path to output file, default to result.data")
parser.add_argument("-b", metavar="NBINS", 
	                        type=int, help="Number of bins")
parser.add_argument("-t", metavar="THRESHOLD", type=float,
	                       help="Threshold, range in (0,1)")
args = parser.parse_args()
print(args)

