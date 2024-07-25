
import pandas as pd
import sys
import shutil
from gt_mat_smartpca.gt_matrix import *

########################################################################

# usage run_sa_qc.py  
raw_geno    = sys.argv[1]
raw_ind     = sys.argv[2] 
filt_log    = raw_geno.split('.geno')[0] + 'filt.log'


########################################################################

n_thres = 0.5

########################################################################
# Write qc log

# Read in names from .ind file
indiv    = read_ind(raw_ind)
df_indiv = pd.DataFrame(indiv)

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
    filt = df.drop(df.columns[bad], axis=1)
    filt_join = filt.apply(lambda x: ''.join(x), axis=1)
    l0 = filt_join.tolist()

if len(bad) == 0:
    filt_ind = raw_ind.split('.ind')[0] + '_filt' + '.ind'
    shutil.copyfile(raw_ind, filt_ind)

# Write filtering log
with open(filt_log, 'w') as f:
    f.write('Filtering on file: ' + raw_geno +'\n')
    f.write('Positions_input: ' + str(df.shape[0]) + '\n')
    f.write('Samples_input: ' + str(df.shape[1]) + '\n')
    f.write('Samples_input: ' + str(df.shape[1]) + '\n')
    f.write('Samples_excluded: ' + str(bad) + '\n')

