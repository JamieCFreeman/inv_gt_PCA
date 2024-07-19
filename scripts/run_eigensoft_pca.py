
import sys

# First write .par file
# p1 = int( sys.argv[1] )
pattern_in  = sys.argv[1]
pattern_out = sys.argv[2]
#pattern_in = '2Lt_mat_filt'
#pattern_out = '2Lt_test_03'
par_file    = pattern_out + '_par.txt'
#par_file    = pattern_out + '/' + pattern_out + '_par.txt'


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

# YAML doesn't accept \t characters!
with open('rmd.temp', 'w') as f:
    f.write('---' + '\n')
    f.write(f"{'title:' :<25}" + '"' + pattern_out + "inversion genotyping by PCA" + '"' + '\n')
    f.write('date: ' + '"`r format(Sys.Date(),\'%B %e, %Y\')`"' + '\n')
    f.write(f"{'output:' :<25}" + '\n')
    f.write( '  ' + f"{'pdf_document:' :<25} default" + '\n')
    f.write(f"{'params:' :<25}" + '\n')
    f.write( '  ' + f"{'genotype_file:' :<25} {pattern_in}.geno" + '\n')
    f.write( '  ' + f"{'snp_file:' :<25} {pattern_in}.snp" + '\n')
    f.write( '  ' + f"{'indiv_file:' :<25} {pattern_in}.ind" + '\n')
    f.write( '  ' + f"{'log_file:' :<25} {pattern_out}.log" + '\n')
    f.write( '  ' + f"{'evec_file_ann:' :<25} {pattern_out}_inv.tsv" + '\n')
    f.write('---' + '\n')

