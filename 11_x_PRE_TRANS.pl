# OUTPUT: each variant in as: 
# 1. reference (in Translation/$C folder) and 
# 2. variant file (in Translation/$C folder/Comparsion folder), chromosome specifically
# The Ref and Alt sequences are visible in lower case letters
#
# ARCHIVED FILES: none
# Elements needed:
# [20] - chromosome
# [21] - variant position
# [33] - Nmer
# [37] - Strand
# [38] - TxStart
# [39] - TxEnd
# [40] - cdsStart
# [41] - cdsEnd
# [42] - Exon count

#!/usr/bin/perl -w
# use 5.020;
# use strict;
use warnings;

use Data::Dumper;
use File::Copy;
use Parallel::ForkManager;
use POSIX;

my$pm=new Parallel::ForkManager(4);

my$T_R='';
my$CAT_R='';
my$Temp_R='';
my$CHR_R='';
my$A_R='';

# make directories for all 24 chromosomes, and format the files into single-line entries:

# Test thoroughly for error messages
opendir(D_t,"$T_R");
	foreach(my$C=1;$C<=24;$C++){
$pm->start and next;	
	mkdir("$T_R/v\_$C");mkdir("$T_R/v\_$C/Comparison");

	open(I,"$T_R/$C");open(O,">>$T_R/f\_$C");select(O);
	while(<I>){s/$/\|/ if/^\d/;chomp if/^\d/;print;}close(I);close(O);
	unlink("$T_R/$C");rename("$T_R/f\_$C","$T_R/$C");

	open(I,"$T_R/$C");while(<I>){my@V=split/ /;s/\|/\n/;
	# Split chromosome file single entry specifically:
	open(O,">>$T_R/v\_$C/$V[0]");select(O);print;print"\n";}close(I);close(O);unlink("$T_R/$C");

# ----------------------------	
	opendir(D_v,"$T_R/v\_$C");
	my@file=grep{/^\d/}readdir D_v;
		foreach my$file(@file){
			# ------ Make the variant position visible on Ref ------------------
			open(I,"$T_R/v\_$C/$file");open(R,">>$T_R/v\_$C/L\_$file");select(R);
				while(<I>){
					if($.==1){my@U=split/ /;					
						# GLOBS are MUST, otherwise the script does not work
						$V_pos=$U[21];					# variant position
						$TxStart=$U[38]; 				# TxStart
						$R=$U[23];						# upper case Reference sequence					
						$A=$U[24];

						$Substr_Pos=($V_pos-$TxStart-1);# position on $.==2 where the lower case substitution is introduced
						$r=lc($R);						# lower case Reference sequence
						$a=lc($A);						# lower case Alt sequence
						$R_L=length($R)}
					
					# ------- print Ref-visible Tx --------------------
					if($.==2){my$Vis_R=substr($_,$Substr_Pos,$R_L,$r)}print;
					}close(I);close(R);

			# ----- Make visibly-mutated Alt: this step uses the GLOBS computed in the previous step (see above)
			open(I,"$T_R/v\_$C/$file");open(A,">>$T_R/v\_$C/Comparison/$file");select(A);	
			while(<I>){if($.==2){my$Vis_A=substr($_,$Substr_Pos,$R_L,$a)}print;}close(I);close(A);
			unlink("$T_R/v\_$C/$file");rename("$T_R/v\_$C/L\_$file","$T_R/v\_$C/$file");
		}
	close(D_v);

$pm->finish;
}
$pm->wait_all_children;
close(D_t);


















