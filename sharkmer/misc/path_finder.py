# Reads the verbose output of sharkmer to find sub kmer paths
# An example line:
#   TCATAAAAATATATGATATATTAAAAAGAA sub_kmer being extended for node 0. There are 1 candidate kmers for extension. Added sub_kmer CATAAAAATATATGATATATTAAAAAGAAA for new node 18 with edge kmer count 8. Path length is 1. There are now 10 unvisited and 11 non-terminal nodes in the graph. 


import sys
import os
import argparse
import re

def main():
	parser = argparse.ArgumentParser(description="Reads the verbose output of sharkmer to find sub kmer paths")
	parser.add_argument("verbose", help="The verbose output of sharkmer")
	parser.add_argument("seed", help="The starting kmer")
	args = parser.parse_args()

	print(f"File: {args.verbose}")
	print(f"Seed: {args.seed}")


	verbose = open(args.verbose)
	path = []
	path.append(args.seed)
	print(f"{args.seed}")

	for line in verbose:
		match = re.search(r'([CGTA]+) sub_kmer .+ Added sub_kmer ([CGTA]+)', line)
		if match:
			print(line)
			parent = match.group(1)
			child = match.group(2)
			if parent in path:
				path.append(child)
				print(f"{child}")

	verbose.close()

if __name__ == "__main__":
	main()
