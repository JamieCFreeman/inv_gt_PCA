
#INVERSIONS = ['1A', '1Be', '2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']
INVERSIONS = ['2LT', '2RNS', '3LP', '3RK', ]

import os
from datetime import date
import gt_mat_smartpca.get_scatter_int as gs

OUTDIR = os.getcwd()
DATE  = date.today()

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
		expand(f"{{inv}}_PCA_run/{{inv}}_{DATE}.pdf", inv=INVERSIONS),
		expand(f"{{inv}}_PCA_run/{{inv}}_mat.{{ext}}", inv=INVERSIONS, ext=EXT)
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
		f"{{inv}}_PCA_run/{{inv}}_mat.{{ext}}"
	run:
		import gt_mat_smartpca.scatter_gather as sg	
		import shutil
		# 
		# For .geno & .snp files, need to merge, but .ind are all the same
		ext = wildcards.ext
		if ext == 'geno' or ext == 'snp':
			sg.merge_big(list(input), output[0])
		elif ext == 'ind':
			shutil.copy2(input[0], output[0])	

rule sa_qc:
	input:
		ind   = f"{{inv}}_PCA_run/{{inv}}_mat.ind",
		geno  = f"{{inv}}_PCA_run/{{inv}}_mat.geno"
	output:
		sa_qc    = f"{{inv}}_PCA_run/{{inv}}_mat_sample_qc.tsv",
		filt_ind = f"{{inv}}_PCA_run/{{inv}}_mat_filt.ind"
	shell:
		"""
		python run_sa_qc.py {input.geno} {input.ind}
		"""

rule filt_mat:
	input:
		geno = f"{{inv}}_PCA_run/{{inv}}_mat.geno",
		snp  = f"{{inv}}_PCA_run/{{inv}}_mat.snp",
		sa_qc = f"{{inv}}_PCA_run/{{inv}}_mat_sample_qc.tsv"
	output:
		tmp  = temp( f"{{inv}}_PCA_run/{{inv}}_mat_filt.temp"),
		geno = f"{{inv}}_PCA_run/{{inv}}_mat_filt.geno",
		snp  = f"{{inv}}_PCA_run/{{inv}}_mat_filt.snp"
	shell:
		"""
		paste {input.snp} {input.geno} |
		awk '$6!="X" {{print $0}}' | awk '$6!="multi" {{print $0}}' | awk '$7 !~ /3/' > {output.tmp}
		awk '{{print $7}}' {output.tmp} > {output.geno}
		cut -f 1-6 {output.tmp} > {output.snp}
		"""

rule write_PCA_par:
	input:
		geno = f"{{inv}}_PCA_run/{{inv}}_mat_filt.geno",
		snp  = f"{{inv}}_PCA_run/{{inv}}_mat_filt.snp",
		ind  = f"{{inv}}_PCA_run/{{inv}}_mat_filt.ind"
	output:
		par = f"{{inv}}_PCA_run/{{inv}}_{DATE}.par"
	shell:
		"""
		python scripts/run_eigensoft_pca.py {workflow.basedir}/{input.geno} {workflow.basedir}/{output.par}
		"""

rule smartpca:
	input:
		par = f"{{inv}}_PCA_run/{{inv}}_{DATE}.par"
	output:
		log  = f"{{inv}}_PCA_run/{{inv}}_{DATE}.log",
		evec = f"{{inv}}_PCA_run/{{inv}}_{DATE}.evec"
	shell:
		"""
		smartpca -p {input.par} > {output.log}
		"""

rule annot_evec:
	input:
		evec = f"{{inv}}_PCA_run/{{inv}}_{DATE}.evec"
	output:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}_inv.tsv"
	shell:
		"""
		python annot_evec.py {input}
		"""

rule rmd_par:
	input:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}_inv.tsv"
	output:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}rmd.temp"
	params:
		date = DATE
	shell:
		"""
		python scripts/write_rmd_par.py {wildcards.inv}_mat_filt \
			{wildcards.inv}_{params.date} {output}
		"""

rule write_rmd:
	input:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}rmd.temp"
	output:
		rmd = f"{{inv}}_PCA_run/{{inv}}_{DATE}.Rmd"
	params:
		rmd = "scripts/plot_PCA.Rmd"
	conda: "envs/rmd.yaml"
	shell:
		"""
		cat {input} <( tail -n +12 {params.rmd}) > {output.rmd}
#		sed -i 's/In.2L.t/{wildcards.inv}/' {output.rmd}
#		Rscript -e "rmarkdown::render(\"{output.rmd}\")"
		"""

rule rmd_report:
	input:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}.Rmd"
	output:
		f"{{inv}}_PCA_run/{{inv}}_{DATE}.pdf"
	conda: "envs/rmd.yaml"
#	script:
#		f"{{inv}}_PCA_run/{{inv}}_{DATE}.Rmd"
	shell:
		r"""
cat <<'EOF' > {rule}.$$.tmp.R

rmarkdown::render("{input}")

EOF
Rscript {rule}.$$.tmp.R
rm {rule}.$$.tmp.R
        """
