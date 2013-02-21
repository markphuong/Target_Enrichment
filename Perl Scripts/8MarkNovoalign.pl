#!/usr/bin/perl

# This script aligns cleaned reads to the 1) final assemblies and 2) assemblies-in-target, respectively, and calculates specificity of each library. It needs a BED file that 
# determines the target region in each of the contigs in the assemblies-in-target. It will finally produce two bam files for each individual library, one is all the reads that
# are mapped to your contigs; another one is reads that are overlapped with your target regions with at least 1bp. Which bam to use for variant calling is really determined by
# the quality of your de novo assembly. The assembly errors increase at the edges of the contigs and I found that this is inevitable, so that is why I more trust reads that are
# overlapped with targets. However, if you believe that your assembled contigs are clean then you may consider to include reads that are mapped to regions outside your target 
# regions as well (without overlapping with the targets). 

# This program requires you to provide a series of options, so please read the following and novoalign mannual to make sure to understand what you what use. 

# Written by Ke Bi.

#use Your::Power;

use warnings;
use strict;

my @TargetAssemblies = <*finalNuclearContig.fa>;
#my $dir = '/home/analysis/Desktop/Mark_analyses/7ConsensusAssemblies/mtDNAtest';
my $insertSize = '160';
my $insertSizeSTDEV = '60';
my $aScore = '90';

foreach my $assembly (@TargetAssemblies) {
	my $read1 = $1 . "_1p_final_renamed.fastq" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	my $read2 = $read1; $read2 =~ s/1p/2p/;
	my $unpaired = $1 . "_u_final_u_combined.fastq" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
#	my $FinalAssembly = $1 . "_contig.fa.final" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	my $BED = $1 . "_BEDNuclear.txt" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;

########### 1. index your reference using novoindex ############

	my $targetNix = $assembly . ".ndx";
#	my $finalNix = $FinalAssembly . ".ndx";

	system ("novoindex ${targetNix} $assembly");
#	system ("novoindex ${finalNix} $FinalAssembly");

########## 2. generate a list of genomic intervals for variants calling ########
#	my $out = $1 . "_region.txt" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;

#	open(IN, "<$BED");
#	open(OUT, '>', $out);
	#my @region = <IN>;
#	while (<IN>) {
#		chomp($_);
#		my @region = split /\s+/, $_;
#		my $gene = $region[0];
#		my $start = $region[1];
#		my $end = $region[2];
#		print OUT "$gene:$start-$end\t";
#	}
#	close(OUT); close(IN); 

########## 3. alignment using novoalign  ###########

	my $pairedOut1 = $1 . "_outPaired1" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	my $unpairedOut1 = $1 . "_outUnpaired1" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
#	my $pairedOut2 = $1 . "_outPaired2" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
#	my $unpairedOut2 = $1 . "_outUnpaired2" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
 
	system("novoalign -d $targetNix -f $read1 $read2  -i PE $insertSize, $insertSizeSTDEV -t $aScore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $pairedOut1");

	system("novoalign -d $targetNix -f $unpaired -t $aScore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $unpairedOut1");

#	system("novoalign -d $finalNix -f $read1 $read2  -i PE $insertSize, $insertSizeSTDEV -t $aScore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $pairedOut2");

#	system("novoalign -d $finalNix -f $unpaired -t $aScore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > $unpairedOut2");

################# 4. remove reads that are not aligned to the refs #####################

	system("grep -v ZS:Z:NM $pairedOut1 > $pairedOut1.sam");
	system("grep -v ZS:Z:NM $unpairedOut1 > $unpairedOut1.sam");
#	system("grep -v ZS:Z:NM $pairedOut2 >  $pairedOut2.sam");
#	system("grep -v ZS:Z:NM $unpairedOut2 > $unpairedOut2.sam");

################# 5. calculate specificity ###################
	my $num1 = 'num1.txt';
	my $num2 = 'num2.txt';
	my $num3 = 'num3.txt';
	my $num4 = 'num4.txt';
	my $num5 = 'num5.txt';
   
	system ("grep HS $pairedOut1.sam | wc -l > $num1");
	system ("grep HS $unpairedOut1.sam | wc -l > $num2");
	system ("grep \@HS1 $read1 | wc -l > $num3");
	system ("grep \@HS1 $read2 | wc -l > $num4");
	system ("grep \@HS1 $unpaired | wc -l > $num5");
    
	open(IN1, '<',$num1);
	open(IN2, '<',$num2);
	open(IN3, '<',$num3);
	open(IN4, '<',$num4);
	open(IN5, '<',$num5);
    
	my @num1 =<IN1>;
	my @num2 =<IN2>;
	my @num3 =<IN3>;
	my @num4 =<IN4>;
	my @num5 =<IN5>;
   
	my @reads_mapped = (@num1, @num2, @num3, @num4, @num5);

	my $specificity = $1 . "_specificity.txt" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
#	my $lib = $1 if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	open (OUT, '>', $specificity);
   
	my $ratio = ($reads_mapped[0]+$reads_mapped[1])/($reads_mapped[2]+$reads_mapped[3]+$reads_mapped[4])*100;
	printf OUT "%.2f", "$ratio";
#	print OUT "%", "\n";
	close(IN1); close(IN2); close(IN3); close(IN4); close(IN5);

################ 6. Samtools is now in play ##################

	my $sorted_in_target_bams = $1 . "_in_target_sorted" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	my $sorted_all_bams = $1 . "_all_sorted" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;    

	system("samtools view -bS $pairedOut1.sam >  $pairedOut1.bam");
	system("samtools view -bS $unpairedOut1.sam > $unpairedOut1.bam");

	system("samtools merge raw.bam $pairedOut1.bam $unpairedOut1.bam ");

	system("samtools sort raw.bam $sorted_all_bams");
    
	my $sorted_all_bams2 = $sorted_all_bams. '.bam';
    
	system("samtools index $sorted_all_bams2");

	my $pileup1 = $1 . ".pileup" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	#my $pileup2 = $1 . "_targetsOnly.pileup" if $assembly =~ m/(index\S+)_finalNuclearContig.fa/;
	
	system("samtools mpileup -f GStarget.txt $sorted_all_bams2 > $pileup1");
	#system("samtools mpileup -f GStarget.txt $sorted_all_bams2 -l $BED > $pileup2");

	#system("samtools view -h -f $assembly $sorted_all_bams2 $out > reads_overlapped_with_target.sam");
    
	#system("samtools view -bS reads_overlapped_with_target.sam > reads_overlapped_with_target.bam");

	#system("samtools sort reads_overlapped_with_target.bam $sorted_in_target_bams");

	#my $sorted_in_target_bams2 = $sorted_in_target_bams . ".bam";

	#system("samtools index $sorted_in_target_bams2");

	system("rm *outPaired* *outUnpaired* raw.bam num*.txt *.ndx"); 

}

#system ("rm $BED");
