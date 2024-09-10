#!/usr/bin/env python

import sys

sys.path.append(r'/home/jamie/FAS1K_utils')
########################################################################

import pandas as pd
import os
from itertools import compress
from gt_mat_smartpca.scatter_gather import *

########################################################################


########################################################################

def get_inv_bk(i, f = '/home/jamie/FAS1K_utils/inv_bk.tsv'):
	'''
	From file with breakpoints f get the breakpoints and arm of
	the specified inversion
    Inv must be in list ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']
	Returned as dict {arm , break1, break2}
	'''
	bk     = pd.read_table(f)
	
	# Should add checkpoint to make sure inv is formatted as expected
	inv_list = ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']

	arm    = bk[bk['inversion']==i]['arm'].values[0]
	break1 = bk[bk['inversion']==i]['break1start'].values[0]
	break2 = bk[bk['inversion']==i]['break2start'].values[0]
	
	return {'arm': arm, 'break1':break1, 'break2': break2}

def match_arm_ref(a, d):
	'''
	For current arm find matching ref fas1k from a dir
	'''
	ref_bool = [ a in x for x in os.listdir(d) ]
	now      = list(compress(os.listdir(d), ref_bool))
	ref      = os.path.join(d, now[0])
	
	return ref

def get_scatter_int(i, s=10000):
	'''
	From inv and chunk size return scatter intervals
	Inv must be in list ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']
    '''
	# Depending on number of files, can alter chunk size so not hogging memory 
	chunk_size = s
	inv_dict = get_inv_bk(i) 
	
	# Get interval list
	parts  = coord_from_chunk(inv_dict['break1'] - 50000, inv_dict['break1'] + 50000, chunk_size)
	parts += coord_from_chunk(inv_dict['break2'] - 50000, inv_dict['break2'] + 50000, chunk_size)
	
	return parts


########################################################################

# If main, Write to file
if __name__ == '__main__':
	inv     = sys.argv[1]
	int_file = inv + '_intervals.txt'
	# Depending on number of files, can alter chunk size so not hogging memory 
	chunk_size = 10000
	
	inv_dict = get_inv_bk(inv)
	ref       = match_arm_ref(inv_dict['arm'], 'ref_fas1k')
	
	intervals = get_scatter_int(inv)
	
	write_int_list(intervals, int_file)
