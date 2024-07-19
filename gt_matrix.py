
#!/usr/bin/env python3

# 2024-06-06 JCF

# Purpose:
# Functions to generate a genotype matrix from a set of fas1k files


###############################################################

import fas1k_utils as f1k
from itertools import compress
import re

###############################################################
# Functions used to generate the matrix

def flatten_list(l):
	return [item for sublist in l for item in sublist]

def allele_count(l):
	# Each entry string is of length sample + 1 for ref
	# Return a count of non-N alleles
	# First, get non-N alleles
	non_n = set([*l]).difference(set('N'))
	# If no heteozygous sites, allele count is just number of non-N alleles
	if (check_het(l, 0) == 1):
		return len(non_n)
	# If het sites, expand them out to their components and merge with homo sites
	elif (check_het(l, 0) == 2):
		# Next, check for heterozygous sites
		het  = check_het(non_n, 1) 
		exp  = flatten_list([ het_dict[x] for x in het ] )
		
		homoz_set = set(('A', 'T', 'C', 'G', 'N'))
		hom = non_n.intersection(homoz_set)
		
		return len( hom.union( set(exp) ))

def all_n(l):
	# Check whether any samples have called sites, 
	# Returns 0 if no called sites
	m = l[1:len(l)]
	return len( set(m).difference(set('N')) )

def check_het(l, verbosity=0):
	'''
	Are any sites heterozygous? Return ploidy of SNP set
	When verbosity is 0, return integer result of ploidy
	When verbosity is 1, return het sites
	'''
	homoz_set = set(('A', 'T', 'C', 'G', 'N'))
	het_set   = set(('Y', 'R', 'W', 'S', 'K', 'M'))
	snp_set   = set(l)
	
	# Does the snp_set intersect with the heterozygous codes?
	if (verbosity == 0):
		if ( len(snp_set.intersection(het_set)) > 0):
			return 2
		elif ( len(snp_set.intersection(het_set)) == 0):
			return 1
	if (verbosity == 1):
		return list( snp_set.intersection(het_set) )

het_dict= {
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"S": ["G", "C"],
	"W": ["A", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"]
	}

def score_invariant(r, x):
	# If site is invariant, there are two possible calls:
	# N
	if (x == 'N'):
		return '9'
	# homozygous ref (0)
	elif (x == r):
		return '0'

def score_biallelic(r, x):
	if   ( x == 'N'):
		return '9'
	# homozygous ref (0)
	elif ( x == r):
		return '0'
	# het (1)
	elif ( check_het(x) == 2):
		return '1'
	# homozy alt (2)
	elif ( check_het(x) == 1):
		return '2'

def gt_site(l):
	# First take each string, and split into indiv char
	v = [*l]
	# Remove ref allele from front
	r = v.pop(0)
	s = ''
	if ( allele_count(l) == 1 ):
		for i in v:
			s = s + score_invariant(r, i)
	if ( allele_count(l) == 2 ):
		for i in v:
			s = s + score_biallelic(r, i)
	return s

def convert_site(l):
	'''
	For a string of alleles
	'''
	# If all sites are n, return a 9 for each sample
	if ( all_n(l) == 0 ):
		s = '9' * ( len(l) -1 )
		return s
	# Now consider completely homozygous sites
	elif ( allele_count(l) <= 2):
		s = gt_site(l)
		return s
	elif ( allele_count(l) >= 3):
		s = '3' * ( len(l) -1 )
		return s


def expand_set_het(s):
	'''
	For a set of nucleotides, expand heterozygous codes into 
	component nucleotides, and return expanded set of nucleotides
	'''
	het_sites = set(het_dict.keys()).intersection(s)
	if len(het_sites) > 0:
		return s.union(het_dict[ list(het_sites)[0] ]) - het_sites
	elif len(het_sites) == 0:
		return s

def get_snp_line(x, chrom, p1, nt_set):
	'''
    Input is a single line from the f1k_zip and its positioning details
    '''
	pos = x + p1
	snp_name = chrom + '_' + str(pos)
	# Ref f1k entry is the first in the string, alt allele will be any other
	#   non-N nt (tri-nt sites will be filtered later- this would arbitrarily choose one)
	ref = nt_set[0]
	alt = expand_set_het( set( nt_set ) ) - { nt_set[0], 'N' }
	
	# Alt above is a set, need to extract nt entry, or
	# If no alt allele, set alt to X as per Eigenstrat spec
	if len(alt) == 1:
		alt = list(alt)[0]
	elif len(alt) == 0:
		alt = 'X'
	elif len(alt) > 1:
		alt = 'multi'
	
	# Position in morgans- set to 0 for unkown
	pos_m = 0
	
	# Columns are snp_name, chr, genetic pos in M (0 if unknown), pos in bp, ref, var 
	#(For monomorphic SNPs, the variant allele can be encoded as X (unknown) )
	snp_line = snp_name + '\t' + str(f1k.arm_to_int(chrom)) + '\t' + str(pos_m) + '\t' + str(pos) + '\t' + ref +'\t' + alt
	
	return snp_line

def get_pop_code(s, k):
	'''
	For a string s, change case to upper, then check against a list 
	of possible known population codes k. 
	'''
	# re.search for each pop code in the provided string, output -> bool, 
	#	then compress list to matches 
	o = list( compress(k, [ bool( re.search(x, s.upper() ) ) for x in k ]) )
	# if len is 1, then pop code is good, should add some error handling for cases where
	#	mult match or no match
	if len(o) == 1:
		return o[0]
	if len(o) == 0:
		return 'missing'

###############################################################

def f1k_zip(sample_list, ref_fas1k, pos_min, pos_max):
	'''
    For given list of fas1k files, get a list over positions from min to max, where
    each entry is a nt from each file in sample_list
    '''
    # Need to iterate over an unknown number of sequences together- would normally use list comprehension here 
	# First create a list of the read in fas1k sequences
	str_list = [ f1k.extract_fas1k_subseq(pos_min, pos_max, ref_fas1k) ] + [ f1k.extract_fas1k_subseq(pos_min, pos_max, x) for x in sample_list ]
	# Then unpack the list, zip them together, and format zip object as list
	zipped    = list(zip(*str_list))
	# Now use list comprehension to join each position into 1 string-
	# 	resulting list has entries equalling the number of sites with each entry
	# 	a string of ref_allele + all sample alleles
	merge   = [ ''.join(x) for x in zipped ]
	
	return merge

def gen_snp_file(merge, chrom, pos_min):
    '''
    Need snp description file
    '''
    #snp_lines = [ get_snp_line(n, chrom, pos_min, merge[n]) for n in range(0,len(merge) - 1) ]
    snp_lines = [ get_snp_line(n, chrom, pos_min, merge[n]) for n in range(0,len(merge))  ]
        
    return snp_lines

def gen_gt_mat(merge):
    # For each position, convert to eigenstrat 0/1/2/9 format string	
    out = [ convert_site(x) for x in merge ]
    return out


###############################################################

# Functions used to filter the matrix

def filt_list(l, f, *args):
	'''
	Takes a list to filter, a function that will return a bool list representing
	the filter action to be taken, and any arguments for the filtering function
	'''
	# If other arguments are provided, use them
	a = [*args]
	if (len(a) > 0):
		b = [ f(x, a[0]) for x in l ]
	else:
		b = [ f(x) for x in l ]
	return list( compress(l, b) )

# If any rows are same across all samples, they are non-informative for PCA

def set_length_string(s):
	# For a string, split into characters and return set length
	return len( set([*s]) )

def apply_set_length_ge1(s):
	# If only one state, uninformative for PCA
	return set_length_string( s.strip()) != 1

def n_count(s):
	# For a string, split into characters and count the 
	n = sum( [ x == '9' for x in [*s] ] )
	return n

def apply_n_thres(s, t):
	# Filter rows with too much missing data (missing in 1/4 of data)
	return n_count(s) < t

def apply_n_thres2(s):
	# If only 1 state + missing, also uninformative
	b = ( set_length_string( s.strip()) == 2 ) & ( n_count(s) > 0 )
	return not b
