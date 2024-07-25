#!/usr/bin/env python3

###############################################################
# import modules
import sys
import os
from itertools import compress
import importlib.util

# import my modules
sys.path.append(r'/home/jamie/FAS1K_utils')

import fas1k_utils as f1k

###############################################################

# Function to import module from abs path
def import_mod_path(MODULE_NAME, MODULE_PATH):
	spec = importlib.util.spec_from_file_location(MODULE_NAME, MODULE_PATH)
	module = importlib.util.module_from_spec(spec)
	sys.modules[MODULE_NAME] = module
	spec.loader.exec_module(module)
	return module

gt_matrix = import_mod_path("gt_matrix", "/home/jamie/inv_gt_PCA/gt_mat_smartpca/gt_matrix.py")
gs        = import_mod_path("get_scatter_int", "/home/jamie/inv_gt_PCA/gt_mat_smartpca/get_scatter_int.py")

from gt_matrix import *

###############################################################

p1       = int( sys.argv[1] )
p2       = int( sys.argv[2] )
out_now  = sys.argv[3]

#out_file = "1A_mat.txt"
#out_now = out_file.split('.txt')[0] + '_' + str(p1) + '_' + str(p2) + '.geno'

inv = out_now.split('_')[0]
arm = gs.get_inv_bk(inv, f = '/home/jamie/FAS1K_utils/inv_bk.tsv')['arm']

# Writing file with append, so want to make sure it doesn't exist when we start
if ( os.path.isfile( out_now ) ):
    sys.exit( "Output file already exists for chunk:" + out_now )

###############################################################

#print( str(p1) + ' ' + str(p2))
all_inv = ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']

###############################################################

ref = "/home/jamie/FAS1K_utils/ref_fas1k/Reference_Chr" + arm + ".fas1k"
d = ["/home/jamie/DGN_compatible/stock_validation/round2/fas1k", "/home/jamie/DGN_compatible/ZI_N/round2/fas1k",
"/home/jamie/DGN_compatible/FR_N/round2/fas1k",
"/home/jamie/Nexus_diploid_fas1k", "/home/jamie/DGN_compatible/FR_N/round2/fas1k", 
"/raid10/backups/genepool/DPGP2plus/wrap1kb/ZI_inbred_diploid", "/raid10/backups/genepool/DPGP2plus/wrap1kb/FR_diploid",
"/home/jamie/dpgp3_sequences"]

to_exclude = "SD136N"

# For a list of directories get list of files
all_files = []
for x in d:
	dl = os.listdir(x)
	all_files += [ os.path.join(x, y) for y in dl ]


def get_fas1k_now(s, c):
	b =  (c in s ) & ( s.endswith(".fas1k") ) & (to_exclude not in s)
	#b =  (c in s ) & ( s.endswith(".fas1k") ) & ( "diploid" in s) & (to_exclude not in s)
	return b

l = filt_list(all_files, get_fas1k_now, arm)
#l = l[0:100]
# Get line name from file path
names = [ f1k.get_name(x) for x in l]


# List of all included populations
pop_list = ['FR', 'ZI', 'EG', 'SD', 'SP', 'EA', 'EF', 'KM', 'CO', 'UG']

# New libraries have different naming patterns, provide list of pop codes 
#   and check library names against known codes- stop execution and return error
#   if pop code is not matched for all samples
if sum([ 'missing' in  get_pop_code(x, pop_list) for x in names ]) > 0:
    missing  = list(compress(names, [ 'missing' in  get_pop_code(x, pop_list) for x in names ]))
    sys.exit('Error: Check provided list of population codes in sample! These samples do not have a match: ' +  ' '.join(missing) )

sa_codes = [ x + '\t' + 'F' + '\t' +  get_pop_code(x, pop_list) for x in names ]

###############################################################

nt_list = f1k_zip(l, ref, p1, p2)
snp_out = gen_snp_file(nt_list, "2L", p1)
out     = gen_gt_mat(nt_list)

#out = gen_gt_mat(l, ref, p1, p2)

with open(out_now, 'a') as f:
	f.write('\n'.join(out) + '\n')

out_snp = out_now.split('.geno')[0] + '.snp'
with open(out_snp, 'w') as f:
    f. write('\n'.join(snp_out) + '\n')

out_indiv = out_now.split('.geno')[0] + '.ind'
with open(out_indiv, 'w') as f:
    f.write('\n'.join(sa_codes) + '\n')

###############################################################


