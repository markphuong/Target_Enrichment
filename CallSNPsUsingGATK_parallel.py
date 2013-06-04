#!/usr/bin/env python

#REQUIRES: GATK, picard-tools

import os
import sys
import argparse
import multiprocessing

def get_args(): #arguments needed to give to this script
	parser = argparse.ArgumentParser(description="run GATK")

	#forces required argument to let it run
	required = parser.add_argument_group("required arguments") 
	required.add_argument("--map", help="textfile with samples to run and what fasta file to match it to", required=True) #A map file with the sample ID and the fasta file it goes to

	return parser.parse_args()

def align(element):

	ID = element[0]
	fastafile = element[1]

	variables = dict(
	reference = 'mitoLoci_' + fastafile + '.fasta', #your reference or reference sequences
	gatk = '/home/analysis/Downloads/GenomeAnalysisTK-2.1-13-g1706365/GenomeAnalysisTK.jar',
	AddOrReplace= '/home/analysis/Downloads/picard-tools-1.79/AddOrReplaceReadGroups.jar',
	ResultsDir = '/home/analysis/Desktop/Mark_analyses/7ConsensusAssemblies/4_24_2013_mtDNAtesting/test/',
	sample = ID,
	bamFile = ID +'_mito_sorted.bam' )

#1. do something [figure out what it does]
#2. sort the bam file
#3. index the bam file
#4. use GATK to call variants + generate likelihood scores
#5. use GATK to get allele depth information in vcf file

	commands = """
	java -jar {AddOrReplace} INPUT= {bamFile} OUTPUT={ResultsDir}{sample}.rg.bam RGID={sample} RGLB=exon_capture RGPL=illumina RGPU=lane3 RGSM={sample}
	samtools sort {ResultsDir}{sample}.rg.bam {ResultsDir}{sample}.rg.sort
	samtools index {ResultsDir}{sample}.rg.sort.bam
	java -Xmx12g -jar {gatk} -T HaplotypeCaller -R {reference} -I {ResultsDir}{sample}.rg.sort.bam -o {ResultsDir}{sample}_part1.vcf
	java -Xmx8g -jar {gatk} -T VariantAnnotator -A DepthPerAlleleBySample -R {reference} -I {ResultsDir}{sample}.rg.sort.bam --variant {ResultsDir}{sample}_part1.vcf -o {ResultsDir}{sample}_mito_annotated.vcf
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd) 

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




