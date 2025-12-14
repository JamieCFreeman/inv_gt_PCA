#!/usr/bin/env python3

# Purpose:
# For a callset produced by this pipeline, combine to get a 
#       table of just inversion calls

########################################################################

import pandas as pd
import sys
import os
from itertools import compress
from functools import reduce

sys.path.append(r'/home/jamie/FAS1K_utils')
from fas1k_utils import get_name

########################################################################

def read_eig(in_file):
	# Read in .evec file, return as list of lists
	# Eigenvector file is space separated- read in by lines
	eigin = open(in_file, 'r')
	lines = eigin.readlines()
	eigin.close()
	
	# Splitting with no argument splits on whitespace
	l = [ x.split() for x in lines ]
	
	return l

def eig_to_df(l):
	# Take the list output from read_eig, make a df with it
	#  rename column to include inversion, and return calls col	
	header_row   = l.pop(0)
	df           = pd.DataFrame( l, columns = header_row )
	df           = df.set_index("Library_ID")
	
	return df

def pull_calls_col(df, inv):
	# subset to calls for unknowns
	df           = df[df['INV_STATUS']=='NA']
	# Rename calls col & subset
	new_col_name = inv
	df           = df.rename({'INV_CALLS':new_col_name}, axis=1)
	
	return df.loc[:, [new_col_name]]

def merge_l_df(df_l):
	# Merge all dataframes using recursive outer merge with functools reduce
	df = reduce(lambda x, y: x.merge(y, how='outer', left_index=True, right_index=True), df_l)
	return df

def gather_calls_to_df(input_files):
	# Read in the eig data and format as df
	# From list of in files, return df of merged calls
	# 1. Get list of dfs with calls from each inv in the set
	df_l = [ pull_calls_col( eig_to_df( read_eig( x ) ), get_name(x, '_')) for x in input_files ]
	out  = merge_l_df(df_l)	
	
	return out

########################################################################

if __name__ == "__main__":
	# Using expand statement in Snakemake to gather all files (space seperated), 
	#	where first element in argv is script path, so exclude
	input_files   = sys.argv[1:]
	df            = gather_calls_to_df(input_files)

	df.to_csv('out.csv', sep= '\t', na_rep='NA')

########################################################################

