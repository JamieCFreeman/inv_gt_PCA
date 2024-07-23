import sys
sys.path.append(r'/home/jamie/FAS1K_utils')

import os
import shutil
from gt_mat_smartpca.gt_matrix import filt_list


# Scattering..

def coord_from_chunk(p1, p2, window):
	'''
	For splitting genomic intervals for scatter-gather operations.
	From 0-based interval start (p1), end (p2), and window size, return
	nested list of intervals (last window truncated to end position).
	No checking for chr end- make sure your ending position is < end of chr. 
	'''
	c = []
	start = p1
	end = p1
	while (end < p2 ):
		end   = start + window
		c.append( [start, end] )
		start = end 
	c[len(c)-1][1] = p2
	
	return c

def write_int_list(c, out):
	'''
	For a nested interval list like that returned from coor_from_chunk,
	print intervals to a tab-separated text file with one interval per line.
	'''
	# Join start and stop coord by tab, then newline between pairs
	s = '\n'.join([ '\t'.join(map(str, x)) for x in c ])
	
	with open(out, 'w') as f:
		f.write( s + '\n')

# Gathering..

def str_match(s, p):
	b =  ( p in s )
	return b

def get_scatter_files(d, p):
	'''
	For a given directory d and a list of patterns p, return files that 
	contain all patterns provided 
	'''
	dir_list = os.listdir(d)
	for x in p:
		b = filt_list(dir_list, str_match, x)
	return b


def merge_big(in_list, out):
    '''
    For a list of inputs, write one at a time and append
    '''
    with open(out,'wb') as wfd:
        for f in in_list:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)


