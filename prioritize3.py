# prioritize GWAS_DATA -b BIN_SIZE -t THRESHOLD -d DESTINATION

import argparse

parser = argparse.ArgumentParser(description='Post-GWAS Prioritization')
parser.add_argument("gwas", metavar='GWAS_DATA_PATH', 
	                          type=argparse.FileType('r'), 
	                          help="Path to GWAS Data")
args = parser.parse_args()
print args

