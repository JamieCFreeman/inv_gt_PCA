# Dependencies:
I was unable to get smartPCA running from the EIGENSOFT bioconda recipe. I had to download and
compile from source, with the llapcke library linked in the make command as below.
git clone https://github.com/DReichLab/EIG.git
make LDLIBS="-llapacke"
make install

I did have to edit the ploteig file to get it to run.
For the ploteig function, I had to edit the perl shebang line from #!/usr/local/bin/perl to #!/usr/bin/perl, 
and comment out line 145 which calls a function on the developer's file system (which otherwise resulted in error
'Can't exec "/home/np29/bin/fixgreen": No such file or directory at ../bin/ploteig line 145, <FF> line 6.'

Make sure smartpca is in your path (export PATH="/home/jamie/EIG/bin:$PATH").

Generates genotype matrix for inversion breakpoint regions, performs PCA, and plots some output/diagnostics.

Changes to do:
Currently the samples are specified in the script run_gt_mat.py, but should move this out (to config file maybe?).
Summary table for all inv
Fix 1Be- cluster on non-African only, then project other populations? or gt just non-Africa
Pop globbing (Aug samples coming up as UG)


Command for 
scp jamie@marula.genetics.wisc.edu:/home/jamie/inv_gt_PCA/res_out.tar .
