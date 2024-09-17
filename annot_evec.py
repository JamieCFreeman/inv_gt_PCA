#!/usr/bin/env python3

# Purpose:
# For an .evec file produced by smartpca, convert to a tsv and 
#	add an annotation column with gold-standard inversion calls

########################################################################

import pandas as pd
import sys
import os
from itertools import compress

sys.path.append(r'/home/jamie/FAS1K_utils')
from fas1k_utils import get_name

########################################################################

def listfullpath(d):
	# If given string, make a list
	out = []
	if isinstance(d, str):
		d = [d]
	for x in d:
		dl = os.listdir(x)
		out += [ os.path.join(x, y) for y in dl ]
	return out


########################################################################

eig_file       = sys.argv[1]
inv_file       = "/home/jamie/FAS1K_utils/DGN_inv_gt_most.tsv"

#inv_file = "/home/jamie/FAS1K_utils/DGN_inv_gt_most.tsv"
syn_dir  = "/home/jamie/dpgp3_sequences/synth_het"

inv = pd.read_table(inv_file)
inv = inv.set_index("Genome")

########################################################################
# Format names in from this file to match breakpoint file formatting
def reformat_inv(x):
	s = x.replace('In(', '').replace(')', '')
	return s.upper()

inv = inv.rename(reformat_inv, axis=1)

now = get_name(eig_file).upper()

########################################################################
# Add synthethic 3LOk het sequences into known

synt       = listfullpath(syn_dir)
files      = list( compress(synt, [ '.fas1k' in x for x in synt ] ))
file_names = [ x.split('/')[len(x.split('/'))-1].split('_Chr3L')[0] for x in files ]

temp = {'SET': 'synt_het', '1A':'ND', '1BE': 'ND', '2LT': 'ND', '2RNS':'ND', \
                '3LOK': 'INV/ST', '3LP': 'ND', '3RK': 'ND', '3RMO': 'ND', '3RP': 'ND', \
                '3RC': 'ND', '2LIC': 'ND', '2RMAL': 'ND'}

tmp = pd.DataFrame([ temp for x in file_names ] , index = file_names)

inv = pd.concat([inv,tmp])


########################################################################
# Read in .evec file and add known inv calls in a new column
# Eigenvector file is space separated- read in by lines
eigin = open(eig_file, 'r')
lines = eigin.readlines()
eigin.close()

# Splitting with no argument splits on whitespace
sep = [ x.split() for x in lines ]
sep.pop(0)

eig = pd.DataFrame(sep,columns=["Library_ID", "E1", "E2", "E3", "E4", "Pop"])
eig = eig.set_index("Library_ID")

# To get join keeping all rows from left
merge = eig.merge(inv[now], left_index = True, right_index = True, how='left')

out_file = eig_file.split('.evec')[0] + '_inv.tsv'
merge.to_csv(out_file, sep='\t', na_rep='NA')


