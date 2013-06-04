#!/usr/bin/env python

#REQUIRES: GATK, picard-tools

import os
import sys

ID = 'index' 
#number = range(51)
#number = number[1:51]
num = (1)

#for num in number:

variables = dict(
reference = 'mitoRef.fa', #your reference or reference sequences
gatk = '/home/analysis/Downloads/GenomeAnalysisTK-2.1-13-g1706365/GenomeAnalysisTK.jar',
AddOrReplace= '/home/analysis/Downloads/picard-tools-1.79/AddOrReplaceReadGroups.jar',
ResultsDir = '/home/analysis/Desktop/Mark_analyses/7ConsensusAssemblies/8make_mtDNA_genomes/2_13_2013_defaultalignmentscoring/4_4_2013_GATK_changeVCFFormat/',
sample = ID + str(num),
bamFile = ID + str(num)+'_mito_sorted.bam' )

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


