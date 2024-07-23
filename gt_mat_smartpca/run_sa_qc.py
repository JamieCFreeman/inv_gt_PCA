
import pandas as pd
import sys
from gt_mat_smartpca.gt_matrix import read_geno, explode_gt

########################################################################

# usage run_sa_qc.py  
raw_geno    = sys.argv[1]
raw_ind     = sys.argv[2] 
filt_log    = raw_geno.split('.geno') + 'filt.log'


########################################################################

n_thres = 0.5

########################################################################


# Read in names from .ind file
indiv    = read_ind(raw_ind)
df_indiv = pd.DataFrame(indiv)

# Read in gts
strip_gt = read_geno(raw)
exploded = explode_gt(strip_gt)

# Set gt to df and get n proportions, het proportions
df = pd.DataFrame(exploded, columns = df_indiv[df_indiv.columns[0]])
n = (df == '9').sum()
prop_n = n / df.shape[0]

het = (df == '1').sum()
prop_het = het / df.shape[0]

qc = pd.DataFrame({'proportion_n': prop_n, 'proportion_het': prop_het})
qc.to_csv(raw.split('.geno')[0] + '_sample_qc.tsv', sep='\t', na_rep='NA')

