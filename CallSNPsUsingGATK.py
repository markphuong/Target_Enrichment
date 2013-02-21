#!/usr/bin/env python

#REQUIRES: GATK, picard-tools

import os
import sys

ID = 'index' 
number = range(51)
number = number[1:51]

for num in number:

	variables = dict(
	reference = 'mitoRef.fa', #your reference or reference sequences
	gatk = '/home/analysis/Downloads/GenomeAnalysisTK-2.1-13-g1706365/GenomeAnalysisTK.jar',
	AddOrReplace= '/home/analysis/Downloads/picard-tools-1.79/AddOrReplaceReadGroups.jar',
	ResultsDir = '/home/analysis/Desktop/Mark_analyses/7ConsensusAssemblies/mtDNAFinalContigs_2_7_2013/2_12_2013_testGATK/',
	sample = ID + str(num),
	bamFile = ID + str(num)+'_mito_sorted.bam' )

	commands = """
	java -jar {AddOrReplace} INPUT= {bamFile} OUTPUT={ResultsDir}{sample}.rg.bam RGID={sample} RGLB=exon_capture RGPL=illumina RGPU=lane3 RGSM={sample}
	samtools sort {ResultsDir}{sample}.rg.bam {ResultsDir}{sample}.rg.sort
	samtools index {ResultsDir}{sample}.rg.sort.bam
	java -Xmx12g -jar {gatk} -T HaplotypeCaller -R {reference} -I {ResultsDir}{sample}.rg.sort.bam -o {ResultsDir}{sample}_mito.vcf
	""".format(**variables)

	cmd_list = commands.split("\n")
	for cmd in cmd_list:
		os.system(cmd)
