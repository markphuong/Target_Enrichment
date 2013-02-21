#####################################################################################
# a script to clean up assemblies, which are redundant and make them less redundant #
# external dependencies: blat (v34), cap3, cd-hit-est 								#
# written by Sonal Singhal, sonal.singhal1 [at] gmail.com, 29 Dec 2011              #
#####################################################################################

use warnings;
use strict;

my $count = "2";

while ($count < 51) {
	my $file = "index" . $count . "_*";
	my $out = "index" . $count . "_contig.fa"; 
	system("cat $file > $out");
	$count++;
}
