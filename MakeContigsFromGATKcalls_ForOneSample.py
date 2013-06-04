#!/usr/bin/env python

import dendropy
import os
import sys
from collections import defaultdict
from pprint import pprint

#Sample ID
ID = 'index'
num = (1)

sample=ID+str(num)

statement = 'working through ' + sample

print statement

VCF = sample + '_mito_annotated.vcf' #how does your VCF file end?

dna = dendropy.DnaCharacterMatrix.get_from_path("mitoRef.fa", "fasta") #Load your reference files!

#Make a dictionary that uses gene name as your key, and the sequence as the value
genes = defaultdict(dict)
for line in dna:
	seq = dna[line]
	genes[str(line)] = str(seq)

#Make a dictionary that uses gene name as your key, and calls upon a list of the bases [in order] of the sequence
refgenes = defaultdict(list)
for gene, seq in genes.iteritems():
	for base in seq:
		refgenes[gene].append(base)

#collect information about variable sites
VCFfile = open(VCF, 'r')
variant = defaultdict(dict)
for line in VCFfile:
	if not line.startswith('#'):
		#get depth
		info = line.split('\t')

		chopinfo = info[7]
		chopinfo = chopinfo.split(';')
		DP = chopinfo[5].split('=')
		#get Genotype and calculate frequency of ALTERNATIVE allele
		GT = info[-1] #grab the last item in list
		GT = GT.split(':') #split up the shit in into a list
		GT = GT[1] #grab the supposed AD scores
		GT = GT.split(',')

		if len(GT) == 2: #this checks if the AD box is actually there, which should always be two numbers! one for the reference allele, one for the alternative allele
			GT = [float(i) for i in GT]
			Freq = GT[1] / (GT[0] + GT[1])
			#makes the following dictionary: {gene => {position => {info about base @ that position}}}		
			variant[info[0]][info[1]] = dict(
				RefBase = info[3],
				AltBase = info[4],
				Depth = DP[1],
				Frequency = Freq
				)

outfile = sample + '_mitoseq.fa'
out = open(outfile, 'w')

# Mask low coverage variants and heterozygous positions with N, and change all correct SNPs

for gene in refgenes: # for gene in your dictionary where the bases are in a list
	if gene in variant: #if the gene has variants
		for position in variant[gene]: #for the positions where there are variants 
			#and they are represented by at least 80% of the reads for that base and have a read depth greater than 5
			if variant[gene][position]['Frequency'] > 0.8 and \
				variant[gene][position]['Depth'] > 5 : 
				 refgenes[gene][int(position)-1] = variant[gene][position]['AltBase'] #change that position to the alternative allele

			elif variant[gene][position]['Frequency'] < 0.2 : 
				continue

			else:
				refgenes[gene][int(position)-1] = 'N' #otherwise, mask it as N

		#then write the sequence to a fasta file
		gene_name = '>' + sample + '_' + gene + '\n'		
		finalseq = ''.join(refgenes[gene]) + '\n'

		out.write(gene_name)
		out.write(finalseq)
		
	else:
		#if there are no variants, just push out the sequence
		gene_name = '>' +sample + '_' + gene + '\n'		
		finalseq = ''.join(refgenes[gene]) + '\n'

		out.write(gene_name)
		out.write(finalseq)

out.close()







