#!/usr/bin/env python3

###############################################################

# import modules
import os
import sys
import subprocess

from gt_mat_smartpca.gt_matrix import filt_list

###############################################################

# Files to gt
#uk = ["/home/jamie/DGN_compatible/stock_validation/round2/fas1k", 
#      "/raid10/jamie/FR_N/round2/fas1k",
#      "/home/jamie/DGN_compatible/ZI_N/round2/fas1k" ]

# Known set
#k = [ "/home/jamie/Nexus_diploid_fas1k",
#      "/raid10/backups/genepool/DPGP2plus/wrap1kb/ZI_inbred_diploid",
#      "/raid10/backups/genepool/DPGP2plus/wrap1kb/FR_diploid",
#      "/home/jamie/dpgp3_sequences",
#      "/home/jamie/dpgp2_sequences",
#      "/raid10/jamie/diploid_fas1k_nomask/CLARK",
#      "/home/jamie/dpgp3_sequences/synth_het"]

# For a list of directories get list of files
def listfullpath(d):
    # If a string is input, it's probably just one dir, make a single entry list
    if isinstance(d, str):
        d = [d]
    out = []
    # For each directory in list, add full path of all files to output list 
    for x in d:
        dl = os.listdir(x)
        out += [ os.path.join(x, y) for y in dl ]
    return out

def sha_return(f):
    '''
    Use subprocess to run sha1sum for a file and return the string
    '''
    return subprocess.check_output(["sha1sum", f ]).decode(sys.stdout.encoding).split(' ')[0]

def match_filetype(s,t):
    '''
    For a file name string, return bool for whether it matches type t.
    '''
    return os.path.splitext(s)[1] == t

def f1k_from_dir(d):
    all_files = listfullpath(d)
    f1k_files = filt_list(all_files, match_filetype, '.fas1k')
    f1k_files.sort()
    return f1k_files

########################################################################

# If main, Write to file
if __name__ == '__main__':
    uk = sys.argv[1].split(',')
    k  = sys.argv[2].split(',')
    
    # Write file list to file so we can check to see if the set 
    #  of .fas1k files in requested directories have changed
    for_hash = f1k_from_dir(uk) + f1k_from_dir(k)
    for_hash.sort()
    with open( 'file_list.txt', 'w') as f:
      f.write('\n'.join(for_hash) + '\n')



