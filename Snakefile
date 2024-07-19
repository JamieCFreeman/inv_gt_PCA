
#INVERSIONS = ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']
INVERSIONS = ['2RNS', '3LP']

import os
import get_scatter_int as gs

OUTDIR = os.getcwd()

def scatter_files(i, e):
        mat_file = i + "_mat" + e
        pat      = mat_file.split(e)[0]
        intervals = gs.get_scatter_int(i)
        scatter_output = [ pat + '_' + str(x[0]) + '_' + str(x[1]) + e for x in intervals ]
        return scatter_output

#all = scatter_files(INVERSIONS[3]) + scatter_files(INVERSIONS[4])

EXT = ['geno', 'ind', 'snp']

rule all:
	input:
		expand(f"{{inv}}_run/{{inv}}.par", inv=INVERSIONS),
		expand(f"{{inv}}_mat.{{ext}}", inv=INVERSIONS, ext=EXT)
#		expand(f"{OUTDIR}/{{inv}}_intervals.txt", inv=INVERSIONS)

rule run_scattered_gt_mat:
	output:
		geno = temp( f"{{inv}}_mat_{{start}}_{{end}}.geno"),
		snp  = temp( f"{{inv}}_mat_{{start}}_{{end}}.snp"),
		ind  = temp( f"{{inv}}_mat_{{start}}_{{end}}.ind")
	shell:
		"""	
		python run_gt_mat.py {wildcards.start} {wildcards.end} {output.geno}
		"""

rule gather_scatter:
	input:
		lambda wildcards: scatter_files(wildcards.inv, '.' + str(wildcards.ext))
	output:
		f"{{inv}}_mat.{{ext}}"
	run:
		import scatter_gather as sg	
		import shutil
		# 
		# For .geno & .snp files, need to merge, but .ind are all the same
		ext = wildcards.ext
		if ext == 'geno' or ext == 'snp':
			sg.merge_big(list(input), output[0])
		elif ext == 'ind':
			shutil.copy2(input[0], output[0])	

rule filt_mat:
	input:
		geno = f"{{inv}}_mat.geno",
		snp  = f"{{inv}}_mat.snp"
	output:
		tmp  = temp( f"{{inv}}_mat_filt.temp"),
		geno = f"{{inv}}_run/{{inv}}_mat_filt.geno",
		snp  = f"{{inv}}_run/{{inv}}_mat_filt.snp"
	shell:
		"""
		paste {input.snp} {input.geno} |
		awk '$6!="X" {{print $0}}' | awk '$6!="multi" {{print $0}}' | awk '$7 !~ /3/' > {output.tmp}
		awk '{{print $7}}' {output.tmp} > {output.geno}
		cut -f 1-6 {output.tmp} > {output.snp}
		"""

rule write_smart_pca:
	input:
		geno = f"{{inv}}_run/{{inv}}_mat_filt.geno",
		snp  = f"{{inv}}_run/{{inv}}_mat_filt.snp"
	output:
		par = f"{{inv}}_run/{{inv}}.par"
	shell:
		"""
		touch {output}
		"""



