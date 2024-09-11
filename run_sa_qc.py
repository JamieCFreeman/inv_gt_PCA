
import pandas as pd
import sys
import shutil
from gt_mat_smartpca.gt_matrix import *

########################################################################

# usage run_sa_qc.py  
raw_geno    = sys.argv[1]
raw_ind     = sys.argv[2] 
filt_log    = raw_geno.split('.geno')[0] + 'filt.log'
filt_ind    = raw_ind.split('.ind')[0] + '_filt' + '.ind'
filt_geno   = raw_geno.split('.geno')[0] + '_filt.tmp'

########################################################################

n_thres = 0.5

########################################################################
# Write qc log

# Read in names from .ind file
indiv    = read_ind(raw_ind)
df_indiv = pd.DataFrame(indiv)
df_indiv = df_indiv.set_index(0, drop=False)

# Read in gts
strip_gt = read_geno(raw_geno)
exploded = explode_gt(strip_gt)

# Set gt to df and get n proportions, het proportions
df = pd.DataFrame(exploded, columns = df_indiv[df_indiv.columns[0]])
n = (df == '9').sum()
prop_n = n / df.shape[0]

het = (df == '1').sum()
prop_het = het / df.shape[0]

qc = pd.DataFrame({'proportion_n': prop_n, 'proportion_het': prop_het})
qc.to_csv(raw_geno.split('.geno')[0] + '_sample_qc.tsv', sep='\t', na_rep='NA')

########################################################################
# Remove samples that don't meet threshold
bad = prop_n.loc[prop_n > n_thres].index.tolist()

if len(bad) > 0:
# Rejoin strings by position to get back to .geno file type
# Write tmp .geno file filtered for individual
    filt = df.drop(labels=bad, axis=1)
    filt_join = filt.apply(lambda x: ''.join(x), axis=1)
    l0 = filt_join.tolist()
    with open(filt_geno, 'w') as f:
         f.write('\n'.join(l0) + '\n')
    # write filtered ind file
    df_indiv  = df_indiv.drop(labels=bad, axis=0)
    filt_join = df_indiv.apply(lambda x: '\t'.join(x), axis=1)
    l0 = filt_join.tolist()
    with open(filt_ind, 'w') as f:
        f.write('\n'.join(l0) + '\n')

# If no samples filtered, just copy original files over
if len(bad) == 0:
    shutil.copyfile(raw_ind, filt_ind)
    shutil.copyfile(raw_geno, filt_geno)

# Write filtering log
with open(filt_log, 'w') as f:
    f.write('Filtering on file: ' + raw_geno +'\n')
    f.write('Positions_input: ' + str(df.shape[0]) + '\n')
    f.write('Samples_input: ' + str(df.shape[1]) + '\n')
    f.write('Samples_input: ' + str(df.shape[1]) + '\n')
    f.write('Samples_excluded: ' + str(bad) + '\n')

