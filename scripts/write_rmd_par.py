
import sys

# First write .par file
# p1 = int( sys.argv[1] )
pattern_in  = sys.argv[1]
pattern_out = sys.argv[2]
out_file    = sys.argv[3]
#pattern_in = '2Lt_mat_filt'
#pattern_out = '2Lt_test_03'
#par_file    = pattern_out + '/' + pattern_out + '_par.txt'



# YAML doesn't accept \t characters!
with open(out_file, 'w') as f:
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

