#!/usr/bin/env python

#SCRIPT TO: pull mtDNA genes from mtDNA genome

import os
import sys
import dendropy

#Sample ID
ID = 'index'
number = range(51)

number = number[1:46]

lociFile = dendropy.DnaCharacterMatrix.get_from_path("mitoLoci_atricapillus.fasta", "fasta") #Load file with reference sequences

#I did this to have a lineage identifier for my samples
lineageDef = {}
for key in [43, 44, 45]:
	lineageDef[key] = 'atricapillus'
for key in [4, 6, 8, 13, 14, 18, 19, 20, 23, 24, 25, 34, 35, 36, 46, 50]:
	lineageDef[key] = 'central'
for key in [42]:
	lineageDef[key] = 'lateralis'
for key in [7, 10, 11, 12, 22, 27, 28, 33, 47, 48, 49]:
	lineageDef[key] = 'northern'
for key in [9, 31, 32]:
	lineageDef[key] = 'variegatus'
for key in [1, 2, 3, 5, 15, 16, 17, 21, 26, 29, 30, 37, 38, 39, 40, 41]:
	lineageDef[key] = 'southern'

#for each locus in your reference
for locus in lociFile:

	print 'Doin\' ' + str(locus)
	outfile = str(locus) + '.fa'
	out = open(outfile, 'w')
	#then loop through each individual and pull out that locus
	for num in number:
		sample = ID + str(num)
		#create a file to store sequences related to that locus [sequences will be pulled from all individuals because
		#sequences now are grouped by individual, not by locus

#		blastOutFile = sample + '_mitoBlastOutput' #which requires the blast file to tell you where the locus is
		seqInFile = sample + '_mitoseq.fa' #and the file to pull the sequence from

		seqFile = dendropy.DnaCharacterMatrix.get_from_path(seqInFile, "fasta") #load up the sequence file

#		blastInfile = open(blastOutFile, 'r') #open the blast file

		for marker in seqFile: #for each line in the blast file
			if str(marker) == str(locus): #if the line in the blast file matches the locus we are currently on
				seq = seqFile[marker]
				seq = str(seq) #get sequence string
#				start = int(info[6]) - 1 #where the sequence starts
#				end = int(info[7]) #where the sequence ends in the mtDNA geome

				outName = '>' + sample + '_' + str(locus) + '_' + lineageDef[num] + '\n' #make a fasta file with the sample index #, the locus, and the lineage identity		
				finalseq = seq + '\n' #past the sequence correctly truncated!

				out.write(outName)
				out.write(finalseq)
			else:
				continue

	out.close()
