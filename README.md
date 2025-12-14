*Drosophila melanogaster* has multiple large, common inversions that segregate at intermediate frequencies.
Commonly, we want to know the inversion status of lines we are planning to use in 
experiments (because the inverisons are recombination modifiers). This pipeline: 
1. Generates a 0/1/2 genotype matrix for the 100 kb region around the breakpoints of each inversion in EIGENSOFT's SMARTPCA compatible format (with .snp and .ind files needed to run EIGENSOFT). This step is scattered, so benefits from giving Snakemake multiple cores to run.
2. Filters the matrix to include only informative sites and remove samples with too much missing data.
3. Runs smartpca by inversion and outputs eigenvalues from the first 10 PCs.
4. Uses set of files of known inversion status to call unknown samples. Currently, calculates distance in the first 4 dimensions, finds the three closest samples of known genotype, and if those are all of one inversion status (std/std, std/inv, or inv/inv), assign the focal sample that genotype. If the three closest samples, vary, call NA.
5. Outputs a .tsv file of inversion calls and a set of R Markdown reports with PC plots and diagnostic information. Check these reports to make sure the clustering of samples looks good before using your callset! 

### Dependencies:
I was unable to get smartPCA running from the EIGENSOFT bioconda recipe. I had to download and
compile from source, with the llapcke library linked in the make command as below.
> git clone https://github.com/DReichLab/EIG.git
make LDLIBS="-llapacke"
make install

I did have to edit the ploteig file to get it to run.
For the ploteig function, I had to edit the perl shebang line from #!/usr/local/bin/perl to #!/usr/bin/perl, 
and comment out line 145 which calls a function on the developer's file system (which otherwise resulted in error
'Can't exec "/home/np29/bin/fixgreen": No such file or directory at ../bin/ploteig line 145, <FF> line 6.'

Make sure smartpca is in your path (export PATH="/home/jamie/EIG/bin:$PATH").

### Changes to do:
-Fix 1Be- cluster on non-African only, then project other populations? or gt just non-Africa
-Would be cool to do something fancier to call gts- k-means clustering? Current method is conservative. Initial tests with clustering seem like it would take a good bit of time, and current method works. 

### To run with Snakemake:
Run with conda:
snakemake --use-conda --conda-frontend conda -c 24

All results are tarred together for easy transfer:
scp jamie@marula.genetics.wisc.edu:/home/jamie/inv_gt_PCA/res_out.tar .
