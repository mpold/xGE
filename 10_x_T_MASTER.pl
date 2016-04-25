# ARCHIVED FILES: none
# Elements needed:
# [32] - Nmer
# [19] - chromosome
# [20] - variant position
# [37] - TxStart
# [38] - TxEnd
# [36] - Strand

# OUTPUT: 2 lines per variant:
# 1. Entire variant description from previous step front-tagged with $. - unique number (line number from C\_$FILE)
#    - total of elements in the the above line now is 34 ([0]-[33] + UCSC GB elements)
#    - Nmer now is [33]
# 2. Entire Tx corresponding to UCSC Tx-coordinates

#!/usr/bin/perl -w
# use strict;
use warnings;

use Data::Dumper;
use File::Copy;
use Parallel::ForkManager;
use POSIX;

# the number in parentheses denotes the number of parallel processes started
# increasing the number of parallel processes much over 4 does not improve performance on 8-core processor
# Rule of thumb: number of parallel processes ~ (number of cores)/2
# number of SNV-translation folders created
my$pm=new Parallel::ForkManager(4);

my$T_R='';
my$CAT_R='';
my$Temp_R='';
my$CHR_R='';
my$A_R='';

# ---------- write parallel scripts to fetch transcripts from genome ----------
opendir (D_c,"$CAT_R");opendir(D_temp,"$Temp_R");
	# Define grep regex as needed:
	my @file=grep{/^(\w)/}readdir D_c;
		foreach my$file(@file){$FILE=$file;open(I,"$CAT_R/$file");
		# make output files and write the top of each chromosome specific script
		foreach(my$i=1;$i<=24;$i++){open(O,">>$Temp_R/$i\.pl");select(O);
			print"\#!\/usr\/bin\/perl \-w\n";
			print"use Data\:\:Dumper;\n";
			print"open(I,\"C:\/ALL_HUMAN_CHROMOSOMES\/chr$i\.fa\");";
			print"open(O,\"\>\>$T_R\/$i\"\);select(O);\nwhile(\<I\>){\n";					
			close(O)}

			while(<I>){
			my@V=split/\t|\,/;
			my$Chr=$V[19];
			my$Tx_length=($V[38]-$V[37]);
			# write the transcript sub-string commands
			foreach(my$C=1;$C<=24;$C++){open(O,">>$Temp_R/$C\.pl");
				if($C eq $Chr){open(O,">>$Temp_R/$C\.pl");select(O);
					print"if(my\$M_$.=substr(\$_,$V[37],$Tx_length)){";
					print"print'$. @V';";	# $. tags each variant with a unique number (line number from C\_$FILE)
					print'print"$';print'M_';print"$.\\n\";";print"\n}\n";
				close(O)}}}close(I);
	foreach(my$j=1;$j<=24;$j++){open(O,">>$Temp_R/$j\.pl");select(O);print "\}close\(I\);close\(O\);\n";}close(O)}
close(D_c);close(D_temp);

# ------------ Execute the parallel perl scripts -----------------------------------------------------------------------
opendir(D_temp,"$Temp_R");
	@file=grep{/\.pl$/}readdir D_temp;
	foreach my$file(@file){
	$pm->start and next;
	
	# delete chromosome specific perl files that contain no variants - prevents perl from opening large chromosome files..
	# ..samples with only few variants will run faster
	open(P,"$Temp_R/$file");$count=0;while(<P>){$count++if/\=substr/}close(P);unlink("$Temp_R/$file")if$count==0;
	system("$Temp_R/$file");
	unlink("$Temp_R/$file");
	unlink("$Temp_R/$file")if$file=~/^\d{1,}\.pl$/;

$pm->finish;
}
$pm->wait_all_children;
close(D_temp);


# -------------- SUMMARIZE FOLDER -----------------------------------------------------------------
opendir(D_t,"$T_R");opendir(D_a,"$A_R");
open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);print"\nChromosome totals before translation\n";
print"Sum of chromosome totals must equals the total of C\_$FILE variants from the previous step\n";
close(O);

@file=grep{/\d/}readdir D_t;
	foreach my $file(@file){open(I,"$T_R/$file");open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
	while(<I>){}my $L_count=$./2;print"$file,$L_count\n";close(I);close(O)}
close(D_t);close(D_a);


























