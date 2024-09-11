

INVERSIONS = ['1A', '1Be', '2LT', '2RNS', '3LP', '3LOK', '3RK', '3RMO', '3RP']
#INVERSIONS = ['2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP' ]

# Some inversions cosmopolitan, some endemic
COSMOP  = ['2LT', '2RNS', '3LP', '3RK', '3RMO', '3RP']
ENDEMIC = ['1A', '1Be', '3LOk']

import os
from datetime import date
import gt_mat_smartpca.get_scatter_int as gs


#######################################################################################
# Want snakemake to trigger rerun when input files used for PCA change (eg add new known
# to improve calling or run new unknown data)
# file_list_hash writes list of files to 'file_list.txt'- hash is on the file list, not
# the files themselves!
import gt_mat_smartpca.file_list_hash as fh

HASH = fh.sha_return('file_list.txt')[0:7]


#######################################################################################


OUTDIR = os.getcwd()
DATE   = date.today()
#DATE = '2024-08-28'

def scatter_files(i, e):
        mat_file = i + "_mat" + e
        pat      = mat_file.split(e)[0]
        intervals = gs.get_scatter_int(i)
        scatter_output = [ pat + '_' + str(x[0]) + '_' + str(x[1]) + e for x in intervals ]
        return scatter_output



EXT = ['geno', 'ind', 'snp']

#######################################################################################

rule all:
	input:
		"res_out.tar",
		expand(f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.pdf", inv=INVERSIONS),
		expand(f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat.{{ext}}", inv=INVERSIONS, ext=EXT)
#		expand(f"{OUTDIR}/{{inv}}_intervals.txt", inv=INVERSIONS)

rule run_scattered_gt_mat:
	output:
		geno = temp( f"{{inv}}_mat_{{start}}_{{end}}.geno"),
		snp  = temp( f"{{inv}}_mat_{{start}}_{{end}}.snp"),
		ind  = temp( f"{{inv}}_mat_{{start}}_{{end}}.ind")
	params: 
		arm =  lambda wildcards: gs.get_inv_bk(wildcards.inv, f = '/home/jamie/FAS1K_utils/inv_bk.tsv')['arm']
	shell:
		"""	
		python run_gt_mat.py {params.arm} {wildcards.start} {wildcards.end} {output.geno}
		"""

rule gather_scatter:
	input:
		lambda wildcards: scatter_files(wildcards.inv, '.' + str(wildcards.ext))
	output:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat.{{ext}}"
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
		ind   = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat.ind",
		geno  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat.geno"
	output:
		sa_qc     = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_sample_qc.tsv",
		filt_ind  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.ind",
		filt_geno = temp(f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.tmp")
	shell:
		"""
		python run_sa_qc.py {input.geno} {input.ind}
		"""

rule filt_mat:
	input:
		geno = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.tmp",
		snp  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat.snp",
		sa_qc = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_sample_qc.tsv"
	output:
		tmp  = temp( f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.temp"),
		geno = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.geno",
		snp  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.snp"
	shell:
		"""
		paste {input.snp} {input.geno} |
		awk '$6!="X" {{print $0}}' | awk '$6!="multi" {{print $0}}' | awk '$7 !~ /3/' > {output.tmp}
		awk '{{print $7}}' {output.tmp} > {output.geno}
		cut -f 1-6 {output.tmp} > {output.snp}
		"""

rule write_PCA_par:
	input:
		geno = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.geno",
		snp  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.snp",
		ind  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_mat_filt.ind"
	output:
		par = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.par"
	shell:
		"""
		python scripts/run_eigensoft_pca.py {workflow.basedir}/{input.geno} {workflow.basedir}/{output.par}
		"""

rule smartpca:
	input:
		par = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.par"
	output:
		log  = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.log",
		evec = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.evec"
	shell:
		"""
		smartpca -p {input.par} > {output.log}
		"""

rule annot_evec:
	input:
		evec = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.evec"
	output:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}_inv.tsv"
	shell:
		"""
		python annot_evec.py {input}
		"""

rule rmd_par:
	input:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}_inv.tsv"
	output:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}rmd.temp"
	params:
		date = DATE
	shell:
		"""
		python scripts/write_rmd_par.py {wildcards.inv}_mat_filt \
			{wildcards.inv}_{params.date} {output}
		"""

rule write_rmd:
	input:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}rmd.temp"
	output:
		rmd = f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.Rmd"
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
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.Rmd"
	output:
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.pdf",
		f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}_inv_CALLS.tsv"
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

rule tar_res:
	input:
		rmd = expand(f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}.pdf", inv=INVERSIONS),
		tsv = expand(f"{HASH}/{{inv}}_PCA_run/{{inv}}_{DATE}_inv_CALLS.tsv", inv=INVERSIONS)
	output:
		"res_out.tar"
	shell:
		"""
		mkdir -p res_out_tmp
		cp {input} ./res_out_tmp
		tar -cvf res_out.tar res_out_tmp
		rm -r res_out_tmp
		"""


