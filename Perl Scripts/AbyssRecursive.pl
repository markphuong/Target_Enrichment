use strict;
use warnings;

my @files = <*1p*>;

my @abysskmer = qw(21 31 41 51 61);

my @cevalue = qw(10 20);

my $dir = '/home/analysis/Desktop/Mark_analyses/5ExtraTrimmedReadsFixed/Assemblies';

###################

foreach my $read1 (@files) {
	my $read2 = $read1;
	$read2 =~ s/1p/2p/;
	my $unpaired = $1 . "_u_final_u_combined.fastq" if $read1 =~ m/(index\S+)_1p_final_renamed.fastq/;
	my $name = $1 if $read1 =~ m/(index\S+)_1p_final_renamed.fastq/;	
	
	foreach my $kmer (@abysskmer) {
		foreach my $cevalue (@cevalue) {
			my $out = $name . "_kmer" . $kmer . "_ce" . $cevalue;
			system ("abyss-pe mpirun='/home/analysis/bin/mpirun --hostfile /home/analysis/Desktop/hostfile' np=12 k=$kmer n=10 E=0 c=$cevalue e=$cevalue in='$read1 $read2' se='$unpaired' name=$dir/$out");
			system("rm $dir/*.adj");
			system("rm $dir/*path*");
			system("rm $dir/*.dot");
			#system("mv $dir/*contigs* $dir/contig");
			#system("mv $dir/*stats $dir/contig");
			#system("mv $dir/*-6.fa $dir/contig");
			system("rm $dir/*.hist");
			#system("rm $dir/*.fa");
			system("rm $dir/*.dist");	
		}
	}
	foreach my $kmer (@abysskmer) {
		my $out = $name . "_kmer" . $kmer;
		system ("abyss-pe mpirun='/home/analysis/bin/mpirun --hostfile /home/analysis/Desktop/hostfile' np=12 k=$kmer n=10 E=0 in='$read1 $read2' se='$unpaired' name=$dir/$out");
		system("rm $dir/*.adj");
		system("rm $dir/*path*");
		system("rm $dir/*.dot");
		#system("mv $dir/*contigs* $dir/contig");
		#system("mv $dir/*stats $dir/contig");
		#system("mv $dir/*-6.fa $dir/contig");
		system("rm $dir/*.hist");
		#system("rm $dir/*.fa");
		system("rm $dir/*.dist");
	}
}
