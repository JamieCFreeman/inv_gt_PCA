#!/usr/bin/env python3

###############################################################

# import modules
import os
import sys
import subprocess

###############################################################

# Files to gt
uk = ["/home/jamie/DGN_compatible/stock_validation/round2/fas1k", "/home/jamie/DGN_compatible/FR_N/round2/fas1k",
                "/home/jamie/DGN_compatible/ZI_N/round2/fas1k" ]

# Known set
k = [ "/home/jamie/Nexus_diploid_fas1k",
        "/raid10/backups/genepool/DPGP2plus/wrap1kb/ZI_inbred_diploid",
        "/raid10/backups/genepool/DPGP2plus/wrap1kb/FR_diploid",
        "/home/jamie/dpgp3_sequences", "/home/jamie/dpgp2_sequences",
        "/raid10/jamie/diploid_fas1k_nomask/CLARK",
        "/home/jamie/dpgp3_sequences/synth_het"]

# For a list of directories get list of files
def listfullpath(d):
    out = []
    for x in d:
        dl = os.listdir(x)
        out += [ os.path.join(x, y) for y in dl ]
    return out

def sha_return(f):
    '''
    Use subprocess to run sha1sum for a file and return the string
    '''
    return subprocess.check_output(["sha1sum", f ]).decode(sys.stdout.encoding).split(' ')[0]


uk_files = listfullpath(uk)
k_files  = listfullpath(k)

# Check to see if files in the requested directories have changed
for_hash = uk_files + k_files
for_hash.sort()
with open( 'file_list.txt', 'w') as f:
   f.write('\n'.join(for_hash) + '\n')



