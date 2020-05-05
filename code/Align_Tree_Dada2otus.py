# Quick script to create the fasta file for the alignments

import argparse
from collections import defaultdict
import os
from ete3 import Tree


parser = argparse.ArgumentParser(description="Script to generate alignment and phylogenetic tree from a tabular list of ASVs")


parser.add_argument("-i", "--input_file", help="input file, one line per OTU", required=True)
parser.add_argument("-o", "--output_prefix", help="output prefix to use", required=True)

args = parser.parse_args()

fasta_file = args.output_prefix + ".fasta"
fasta_output = open(fasta_file, 'w')

conversion_table = defaultdict(str)
fasta_header = 1


for line in open(args.input_file, 'r'):
    line = line.rstrip()
    fasta_output.write(">" + str(fasta_header) + "\n")
    fasta_output.write(line + "\n")

    conversion_table[str(fasta_header)] = line
    fasta_header += 1

fasta_output.close()

# Run MAFFT
print("Running alignment with " + fasta_file)
aln_file = args.output_prefix + ".aln"
mafft_cmd = "mafft --auto --thread 20 " + fasta_file + " > " + aln_file
os.system(mafft_cmd)


# Run Fasttree
tree_file = args.output_prefix + ".tre"
print("Running tree with " + tree_file)
fastree_cmd = "FastTree -gtr -nt < " + aln_file + " > " + tree_file
os.system(fastree_cmd)


# Then modify the tree node names

t = Tree(tree_file)

for leaf in t:
    original_name = conversion_table[leaf.name]

    leaf.name = original_name

t.write(outfile=args.output_prefix + ".final.tre")

