
import sys

# First write .par file
# p1 = int( sys.argv[1] )
geno_in   = sys.argv[1]
par_file  = sys.argv[2]
#pattern_in = '2Lt_mat_filt'
#pattern_out = '2Lt_test_03'
#par_file    = pattern_out + '_par.txt'
#par_file    = pattern_out + '/' + pattern_out + '_par.txt'

pattern_in   = geno_in.split('.')[0]
pattern_out  = par_file.split('.')[0]

# User var
# Use alternative normalization style? Default is yes
alt_norm = 'NO'

# Number eigenvectors out (default 2)
eig_out = 4

# Number outlier iterations (0 turns off outlier removal)
outlier_it = 0

grm_out = 'grmjunk'

with open(par_file, 'w') as f:
    f.write(f"{'genotypename:' :<25} {pattern_in}.geno" + '\n')
    f.write(f"{'snpname:' :<25} {pattern_in}.snp" + '\n')
    f.write(f"{'indivname:' :<25} {pattern_in}.ind" + '\n')
    f.write(f"{'evecoutname:' :<25} {pattern_out}.evec" + '\n')
    f.write(f"{'evaloutname:' :<25} {pattern_out}.eval" + '\n')
    f.write(f"{'altnormstyle:' :<25} {alt_norm}" + '\n')
    f.write(f"{'numoutevec:' :<25} {eig_out}" + '\n')
    f.write(f"{'numoutlieriter:' :<25} {outlier_it}" + '\n')
    f.write(f"{'familynames:' :<25} NO" + '\n')
    f.write(f"{'grmoutname:' :<25} {grm_out}" + '\n')


command = "smartpca";
command += " -p " + par_file +  " > " + pattern_out + ".log"
print(command)

