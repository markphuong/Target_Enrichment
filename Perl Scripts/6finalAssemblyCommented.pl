#####################################################################################
# a script to clean up assemblies, which are redundant and make them less redundant #
# external dependencies: blat (v34), cap3, cd-hit-est 								#
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 29 Dec 2011              #
#####################################################################################

use warnings;
use strict;

#directory that contains all assemblies
my $dir = '/home/analysis/Desktop/Mark_analyses/6ContigsToMerge/*_contig.fa';
#get rid of any contigs that are shorter than this 
my $minLength = '200';
#how much memory you can dedicate to these external dependencies, in MB
my $mem = '4000';
#the difference level at which we can start to try to cluster contigs
my $dist = '0.99';
#two contigs need to overlap over $overCut% of their potential overlap in order for us to even try to assemble them
my $overCut = '0.95';

###########################
# run the subroutines     #
###########################

my @assembly = <$dir*>; #put all assemblys from a certain directory into an array
foreach my $assembly (@assembly) { #for each assembly in the array.. 
	my $orig = $assembly . ".original"; #make a file and add ".original" to the end of it
	my $call1 = system("cp $assembly $orig"); #copies original assembly to the $assembly.original file

	#first need to rename contigs, so that all the headers are >contig# (or >contig1, contig2, contig3)
	renameContigs($assembly);
	
	#then need to cluster at 100%; also remove contigs shorter than $minLength, this is now reading the file that has been renamed, the $assembly variable has changed
	cluster100($assembly);
	
	#then cluster at 99%, and use cluster info to run cap3
	my $cl99a = cluster99($assembly);
	
	#cl99a is a new fasta file with nearly identical contigs assembled with cap3
	#then take clustered contigs and do blat 1x at 99%
	my $cl99b = clusterBlat($cl99a,'.cl99b');
	
	$dist = $dist - 0.01;
	#then take clustered contigs and do blat 2x at 98%
	my $cl98a = clusterBlat($cl99b,'.cl98a');
	my $cl98b = clusterBlat($cl98a,'.cl98b');
	
	my $final = $assembly . ".final";
	my $call2 = system("cp $cl98b $final");
	}
		
###########################
# behold the subroutines  #
###########################	
	
sub clusterBlat {
    my ($seq,$index) = @_;
    my $clusters = 'blatSelCheck' . '.out' ;
	my $call = system("blat $seq $seq $clusters -noHead -dots=100");
	
	#make a hash in which the sequence is linked to it's contig#
	my (%seq, $id);
	open(SEQ, "<$seq");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) {
			$id = $1;
			}
		else { 
			$seq{$id} .= $line;
			}
		}
	close(SEQ);
	
	my (%clusters, %revClusters);
	my $tracker = 1;
	open(IN, "<$clusters"); #while reading the blat output file
	while(<IN>){
		chomp(my $line = $_);
		my @d = split(/\t/,$line);
		unless($d[9] eq $d[13]) {
			my $overlap;
			#calculate potential overlap
			if ($d[11] < 0.05*$d[10]) {
				#query is hanging off the right end of the hit
				 if ($d[15] + $d[10] > $d[14]) {
				 	$overlap =  $d[14] - $d[15];
				 	}
				#query is entirely within the hit
				 else {
				 	$overlap = $d[10];
				 	}		
				}
			else {
				#hit is entirely within the query
				if ($d[11] + $d[15] < $d[10]) {
					$overlap = $d[14];
					}
				#query is hanging off the left end of the hit
				else {
					$overlap = $d[10] - $d[11];
					}
				}	
			if (($d[0]+$d[1])/$overlap >= $overCut && $d[0]/($d[0] + $d[1]) >= $dist) {
				#overlap is good enough that I am willing to try to assemble these...

				#complicated note-taking because of weird data-structure i have created...
				my ($query, $hit) = 0;
				$query = $revClusters{$d[9]} if $revClusters{$d[9]};
				$hit = $revClusters{$d[13]} if $revClusters{$d[13]};

				if ($query && $hit) {
					if ($query ne $hit) {
						#need to combine these two clusters
						foreach my $hit_id (keys %{$clusters{$hit}}) {
							$clusters{$query}{$hit_id}++;
							$revClusters{$hit_id} = $query;
							}
						delete($clusters{$hit});
						}
					}
				elsif ($query) {
					$clusters{$query}{$d[13]}++;
					$revClusters{$d[13]} = $query;
					}
				elsif ($hit) {
					$clusters{$hit}{$d[9]}++;
					$revClusters{$d[9]} = $hit;
					}
				else {
					$clusters{$tracker}{$d[9]}++; $clusters{$tracker}{$d[13]}++;
					$revClusters{$d[9]} = $tracker; $revClusters{$d[13]} = $tracker;	
					$tracker++;
					}
				}			
			}
		}
	unlink($clusters);

	#make cluster files here 
	my $out = $1 . $index if $seq =~ m/(.*)\.[a-z|0-9]+/i;
	$tracker = 1;
	open(OUT, ">$out");
	my $temp = "localAssembly.fa"; 
	foreach my $cluster (keys %clusters) {
		open(TEMP, ">$temp");
		foreach my $seqid (keys %{$clusters{$cluster}}) {
			print TEMP ">", $seqid, "\n", $seq{$seqid}, "\n";
			delete($seq{$seqid});
			}
		close(TEMP);
		
		my $call = system("cap3 " . $temp . " -z 1 -o 16 -e 11 -s 251");
		my $assembled = $temp . ".cap.contigs";
		my $singlets = $temp . ".cap.singlets";
		
		open(SIN, "<$singlets");
		while(<SIN>) {
			chomp(my $line = $_);
			if ($line =~ m/>/) {
				print OUT ">contig", $tracker, "\n";
				$tracker++;
				}
			else {
				print OUT $line, "\n";
				}
			}
		close(SIN);
	
		
		open(CON, "<$assembled");
		while(<CON>) {
			chomp(my $line = $_);
			if ($line =~ m/>/) {
				print OUT ">contig", $tracker, "\n";
				$tracker++;
				}
			else { 
				print OUT $line, "\n";
				}
			}
		close(CON);
		}

	foreach (keys %seq) {
		print OUT ">contig", $tracker, "\n", $seq{$_}, "\n";
		$tracker++;
		}
	close(OUT);

	my $call3 = system("rm $temp" . "*");	
	return($out);	
	}
	
sub cluster99 {
	my ($assembly) = @_;
	
	my $out = $assembly . '2'; #make a new outfile
	my $call1 = system("cd-hit-est -i $assembly -o $out -c 0.99 -M $mem -r 1 -d 30 -B 1"); #merge contigs that are 99% similar -d length of description in .clsterfile, so 30 characters I am guessing
	my $clusters = $out . ".clstr"; #make this clster file into a variable
	
	#get all the cluster info; this while function will make the cluster number an ID and assign to those clusters the contig numbers associated in that cluster
	my %clusters; my $c; #initialize these two things
	open(IN, "<$clusters"); #open file handle of the cluster file
	while(<IN>) { #while reading the file IN
		chomp(my $line = $_); #remove new line character from each linhe
		if ($line =~ m/>(Clus.*)/) { #for each line that you do that, if it equals >clustersomething, then, make c = to Cluster #, then change it to Cluster#
			$c = $1; $c =~ s/\s//g; #\s is a whitespace character
			}
		else { #for all other lines you find:
			push(@{$clusters{$c}},$1) if $line =~ m/>([A-Z|0-9]+)/i; # $clusters{$c}, $clusters refers to the has, %clusters, {$c} refers to the key, and $clusters{$c} acts as the reference to an array when @{} is placed around it!! $1 gets added to that array. Perl automatically creates the reference to a new empty array. 
			}	
		}
	close(IN);
	my $call2 = system("rm $out" . "*"); 
	
	#get the sequence info; in the hash %seq, associate the sequence to the key which is coded as the Contig #
	my (%seq, $id);
	open(SEQ, "<$assembly");
	while(<SEQ>) {
		chomp(my $line = $_);
		if ($line =~ m/^>(\S+)/) { #if the line is the header
			$id = $1; #make the variable $id the header
			}
		else { 
			$seq{$id} .= $line; #then in the hash, for each header, associate the ID with it
			}
		}
	close(SEQ);
	
	my $cl99 = $assembly . ".cl99a";
	my $tracker = 1;
	open(OUT, ">$cl99");
	my $temp = "localAssembly.fa"; 
	#now go through each cluster and do a "local assembly"
	foreach my $c (keys %clusters) { #for each key in the %cluster
		my @contigs = @{$clusters{$c}}; #make this @contigs array equal to the referenced array in %clusters
		#loner
		if (scalar(@contigs) == 1) { #counts the number of elements in the array
			print OUT ">contig", $tracker, "\n", $seq{$contigs[0]}, "\n"; # $contigs[0] calls upon the first element in @contigs [and there is only one], then $seq{$contigs[0]} calls upon that contig from the hash seq. prints a new header >contig# - new line - then the sequence
			$tracker++;
			}
		#more than one	
		else {	
			open(TMP, ">$temp"); #open a temporary file to write onto
			foreach my $id (@contigs) { #for each ID in the contigs array (which is now >1), 
				print TMP ">", $id, "\n", $seq{$id}, "\n"; #put them in this temporary file [and they are linked by the Contig ID] though, you rename them in the outfile
				}
			close(TMP);
			
			#now need to assemble this...
			my $call = system("cap3 " . $temp . " -z 1 -o 16 -e 11 -s 251"); # z is minimum number of good reads at clip position -o is the overlap length cutoff, -e is the clearance between no. of different, -s is overlap similarity score cutoff
			my $assembled = $temp . ".cap.contigs"; #declare  
			my $singlets = $temp . ".cap.singlets"; #declare
			
			#I am assuming contigs get merged, some get cut off, and this is a code to write into a new file a new contig with which was not assembled in the previous function
			open(SIN, "<$singlets"); #open this file
			while(<SIN>) {
				chomp(my $line = $_);
				if ($line =~ m/>/) { 
					print OUT ">contig", $tracker, "\n";
					$tracker++;
					}
				else {
					print OUT $line, "\n";
					}
				}
			close(SIN);
			#and this one prints a newly assembled contig
			open(CON, "<$assembled");
			while(<CON>) {
				chomp(my $line = $_);
				if ($line =~ m/>/) {
					print OUT ">contig", $tracker, "\n";
					$tracker++;
					}
				else {
					print OUT $line, "\n";
					}
				}
			close(CON);			
			}	
		}	
	close(OUT);
	my $call3 = system("rm $temp" . "*");	
	return($cl99);
	}

sub cluster100 {
	my ($assembly) = @_;
	
	my $out = $assembly . '2'; #following logic of the code..I should get a file that is $assembly.2.2
	my $call1 = system("cd-hit-est -i $assembly -o $out -c 1.00 -l $minLength -M $mem -r 1 -B 1"); #-i infile -o outfile -c sequence_identity_threshold -l throw_Away_length -M how_much_mem_to_use -r both_directions -B store_work_in_disk_space. This function clusters identical sequences..
	rename($out,$assembly); #deletes the $out file and stores that information in the original $assembly file
	my $call2 = system("rm $out" . "*"); #remove extraneous files created by cd-hit-est
	}	
	
sub renameContigs {
	my ($assembly) = @_; #load in assemblies

	my $out = $assembly . '2'; #make an out file with a "2"
	open(IN, "<$assembly"); #open file handle for the assembly
	open(OUT, ">$out"); #open the $assembly.2 file
	my $contig = 1; #variable to start naming the cotings
	
	while(<IN>) { #while in the file
		chomp(my $line = $_); #remove the new line symbol
		if ($line =~ m/>/) { #if the line has >
			print OUT ">contig", $contig, "\n"; #print ">contig" + the number + a new line
			$contig++; #then move on
			}
		else { #and if contig does not have the >
			print OUT $line, "\n"; #print the line [your sequence]
			}
		}
	close(IN); close(OUT); #close both filehandles
	rename($out,$assembly); #deletes the $out file and stores that information in the original $assembly file
	}	
