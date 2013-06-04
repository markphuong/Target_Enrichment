#!/usr/bin/env python

#This script will create a new fasta file that will have contigs that are modified [from original reference] based on reads of a sample

import dendropy
import os
import sys
from collections import defaultdict
from pprint import pprint
import argparse
import multiprocessing

def get_args(): #arguments needed to give to this script
	parser = argparse.ArgumentParser(description="call snps and make contigs")

	#forces required argument to let it run
	required = parser.add_argument_group("required arguments") 
	required.add_argument("--map", help="textfile with samples to run and what fasta file to match it to", required=True) #A map file with the sample ID and the fasta file it goes to

	return parser.parse_args()

def align(element):

	ID = element[0]
	fastafile = element[1]

	sample=ID
	statement = 'working through ' + sample
	print statement
	reference = 'mitoLoci_' + fastafile + '.fasta' #your reference or reference sequences
	dna = dendropy.DnaCharacterMatrix.get_from_path(reference, "fasta") #Load your reference files!
	VCF = sample + '_mito_annotated.vcf' #how does your VCF file end?

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
			call = GT[0]
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
					Frequency = Freq,
					Genotype = call
					)

	outfile = sample + '_mitoseq.fa'
	out = open(outfile, 'w')

	# Mask low coverage variants and heterozygous positions with N, and change all correct SNPs

	for gene in refgenes: # for gene in your dictionary where the bases are in a list
		try:
			variant
		except NameError:
			gene_name = '>' + gene + '\n'		

			finalseq = ''.join(refgenes[gene]) + '\n'

			out.write(gene_name)
			out.write(finalseq)

		if gene in variant: #if the gene has variants
			for position in variant[gene]: #for the positions where there are variants 
				#and they are fixed for the alternative allele and have a read depth greater than 5
				if variant[gene][position]['Genotype'] == '1/1' and \
					variant[gene][position]['Depth'] > 5 :

					#IF It is an InDel and the reference is LONGER than the alternative
					if len(variant[gene][position]['RefBase']) > len(variant[gene][position]['AltBase']):
						length = len(variant[gene][position]['RefBase'])
						
						#delete the number of bases AFTER the position
						for i in range(length): #range goes from 0 to n
							refgenes[gene][int(position)-1+(i)] = '' #and put nothing in there
						
					refgenes[gene][int(position)-1] = variant[gene][position]['AltBase'] #change that position to the alternative allele

				#or they are represented by at least 80% of the reads for that base and have a read depth greater than 5
				elif variant[gene][position]['Frequency'] > 0.8 and \
					variant[gene][position]['Depth'] > 5 :
					#IF It is an InDel and the reference is LONGER than the alternative
					if len(variant[gene][position]['RefBase']) > len(variant[gene][position]['AltBase']):
						length = len(variant[gene][position]['RefBase'])
						
						#delete the number of bases AFTER the position
						for i in range(length): #minus one because you don't delete the first spot
							refgenes[gene][int(position)-1+(i)] = '' #and put nothing in there

					refgenes[gene][int(position)-1] = variant[gene][position]['AltBase'] #change that position to the alternative allele

				elif variant[gene][position]['Frequency'] < 0.2 : 
					continue
				else:
					refgenes[gene][int(position)-1] = 'N' #otherwise, mask it as N

			#then write the sequence to a fasta file
			gene_name = '>' + gene + '\n'		
			finalseq = ''.join(refgenes[gene]) + '\n'

			out.write(gene_name)
			out.write(finalseq)
		
		else:
			#if there are no variants, just push out the sequence
			gene_name = '>' + gene + '\n'		
			finalseq = ''.join(refgenes[gene]) + '\n'

			out.write(gene_name)
			out.write(finalseq)

	out.close()

def main():
	args = get_args() 

	#Make a list of lists, each list within the list will have the first and second elements of the map file that are separated by a tab
	mylist = []
	with open(args.map) as rfile:
		for line in rfile:
			line = line.strip()
			mylist.append(line.split("\t"))

	#start the multiprocessing
	pool = multiprocessing.Pool()
	pool.map(align, mylist)#run the function with the arguments

if __name__ == "__main__": #run main over multiple processors
	main()






