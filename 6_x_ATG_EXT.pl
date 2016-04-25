# OUTPUT FILES:
# ext_atg_m_BIG - MINUS strand ATG-affectors, the variants with possible upstream ATG extension tagged with the number of triplets and DNA extension sequence
# ext_atg_p_BIG - PLUS    "      " ....
# no_atg_BIG - variants with start ATG not affected
# 
# ARCHIVED FILES:
# 1. START_$FILE - merger of the above three files
# 
# OUTPUT ELEMENTS: 21 + UCSC exome descriptions
# Elements from this script:
# [0] - ATG_EXT - if consensus start ATG is affected then either 'ext_atg_m_BIG' or 'ext_atg_p_BIG' or 'no_atg_BIG'
# [1] - EXT_SEQ - if and upstream in-frame ATG without in-frame stop codon (5' ATG(infr)..no_STOP(infr)..ATG(consensus) then value = DNA seq); if foregoing does not apply then value = ''
# Elements from the previous script
# [2] - previous step: [0] - TAG_1: transcript mapping file prefix
# [3] - previous step: [1] - MAP: computed transcript mapping values as shown above (on the right side)
# [4] - previous step: [2] - AP_SM: distance from either plus strand ATG or minus strand STOP
# [5] - previous step: [3] - AM_SP: distance from either minus strand ATG or plus strand STOP
# [6] - [20] - previous step: [4] - [18] - the [0] - [14] from x_SORT.pl; all elements now tab-separated followed by..
# .. NR- or NM-record from UCSC Genome Browser download


#!/usr/bin/perl -w
#use strict;

use Data::Dumper;
use File::Copy;
use POSIX;
use List::Util qw[min max];


my$I_R='';
my$Temp_R='';
my$A_R='';
my$COM_R='';
my$M_R='';
my$CAT_R='';

opendir(D_m,"$M_R");opendir(D_c,"$CAT_R");
	my@file=grep{/^(\w|\d)/}readdir D_m;
		foreach my$file(@file){
# make output files:		
my @OUT_F=('ext\_atg\_p','ext\_atg\_m','no\_atg');
foreach my$OUT_F(@OUT_F){open(OUT_F,">>$CAT_R/$OUT_F\_$file");close(OUT_F)}
$FILE=$file;
# STEP1: map all variants from that qualify from the previous steps: 
			open(IN,"$M_R/$file");
				while(<IN>){
# no warnings;
					my@U=split/\t|\,/; # first comma denotes the start of UCSC GB data
# PART 1: map out the NM and NR EXON/INTRON placements of variants
# add 11 to all UCSC GB elements
# U[6] - variant
# U[23] - TxStart
# U[24] - TxEnd
# U[25] - cdsStart
# U[26] - cdsEnd
# U[18] - 25mer
# U[22] - strand
# U[2] - variation start correction
# U[17] - 1- or 0-based code (see previous scripts)

# PLUS STRAND ACCEPTOR JUNCTION:
# Excel control formula: =MID(U1,(98 - C1),2)
# Excel file: C:\xGeSOFT\xCOMMANDS\1_4_UNIVERSAL_TRANSLATE\INTERESTING_EXCEPTION\CONTROLS\Variation_Type_Controls_01262016\atg_pC_BIG.xlsx
my$cds_P=(100-$U[2]);
my$ATG_P=substr($U[18],$cds_P, 3);
my$STP_P='(TAG|TGA|TAA)';			# Defines stop codons

						if($U[22]=~/\+/&&/^atg\_/){
						open(O,">>$CAT_R/ext\_atg\_p\_$file");select(O);
							foreach(my$up=$cds_P;$up>=($cds_P-97);$up--){
								my $IN_FR=substr($U[18],$up,3);
								my $UP_ATG=($cds_P-$up)/3;	# defines number of codons for the upstream ATG
								my $UP_Seq=($cds_P-$up);	# defines the number of bases of upstream sequence
								my $EXT=substr($U[18],($cds_P-$UP_Seq),$UP_Seq);	# sequence of interest

								my@N=(1..33);			# defines integers -> if $UP_ATG is not integer then out-of-frame ATG									
								foreach my$N(@N){							
my$INFRAME_STOP_P= 
( # Eliminate extension with inframe STPOP between the upstream ATG and original ATG; OOGGLIE!!!!
($EXT=~/$IN_FR\w{0}$STP_P/)||
($EXT=~/$IN_FR\w{3}$STP_P/)||
($EXT=~/$IN_FR\w{6}$STP_P/)||
($EXT=~/$IN_FR\w{9}$STP_P/)||
($EXT=~/$IN_FR\w{12}$STP_P/)||
($EXT=~/$IN_FR\w{15}$STP_P/)||
($EXT=~/$IN_FR\w{18}$STP_P/)||
($EXT=~/$IN_FR\w{21}$STP_P/)||
($EXT=~/$IN_FR\w{24}$STP_P/)||
($EXT=~/$IN_FR\w{27}$STP_P/)||
($EXT=~/$IN_FR\w{30}$STP_P/)
);
										# upstream inframe start, not followd by inframe stop before the original start
										if(($IN_FR eq'ATG')&&($UP_ATG==$N)&&!$INFRAME_STOP_P){
										print"$UP_ATG\-$EXT";
										}}}print"\t$_"}close(O);
# MINUS STRAND:
# Excel control formula: =MID(S2,98-C2,3) -> 98 in Excel == 97 in perl
my$cds_M=(97-$U[3]);
my$TAC_M=substr($U[18],$cds_M,3);
my$STP_M='(CTA|TCA|TTA)';			# Defines stop codons

						if($U[22]=~/\-/&&/^atg\_/){open(O,">>$CAT_R/ext\_atg\_m\_$file");select(O);
							foreach(my$down=$cds_M;$down<=($cds_M+100);$down++){
								my$in_fr=substr($U[18],$down,3);
								my$DOWN_TAC=($down-$cds_M)/3;				# defines number of codons for the downstream ATG
								my$DOWN_Seq=($down-$cds_M);					# defines the number of bases of downstream sequence
								my$ext=substr($U[18],($cds_M+3),$DOWN_Seq);	# inframe piece of DNA downstream of -strand CAT

								my@n =(1..33);	# defines integers -> if $DOWN_TAC is not integer then out-of-frame ATG
								foreach my$n(@n){
my$INFRAME_STOP_M = 
( # Eliminate extension with inframe STPOP between the downstream CAT and original CAT
($ext=~/$STP_M\w{0}$in_fr/)||
($ext=~/$STP_M\w{3}$in_fr/)||
($ext=~/$STP_M\w{6}$in_fr/)||
($ext=~/$STP_M\w{9}$in_fr/)||
($ext=~/$STP_M\w{12}$in_fr/)||
($ext=~/$STP_M\w{15}$in_fr/)||
($ext=~/$STP_M\w{18}$in_fr/)||
($ext=~/$STP_M\w{21}$in_fr/)||
($ext=~/$STP_M\w{24}$in_fr/)||
($ext=~/$STP_M\w{27}$in_fr/)||
($ext=~/$STP_M\w{30}$in_fr/)
);
									# upstream infraem start, not followd by inframe stop before the original start
									if(($in_fr eq 'CAT')&&($DOWN_TAC==$n)&&!$INFRAME_STOP_M){
									print"$DOWN_TAC\-$ext"}}}print"\t$_"}close(O);

						# ---------------- variants of no ATG affection -------------------------------------------
						if(!/^atg\_/){open(OTHER,">>$CAT_R/no\_atg\_$file");select(OTHER);print"\t$_";close(OTHER)}
# Same start and end rules apply to PLUS and MINUS strand Tx, CDS, and Exon starts and ends
# Intron start is right next to exon end, and Intron end right next to exon start
#  PLUS
#  ------->
#  UTR5->- - - - - - cds - - - - - ->- -UTR3->
#  - - - A T G . . . . . . . . T G A - - - - - 
# -2-1 0 1 2 3 4 5 6 . . . . .-2-1 0 1 2 3 4 5
# 
#  									     MINUS
#  									  <-------
#  - - - T C A . . . . . . . . C A T - - - - - 
# -2-1 0 1 2 3 4 5 6 . . . . .-2-1 0 1 2 3 4 5
# <-UTR3<- - - - - - cds - - - - - -<- - -UTR5
			}close(IN);unlink("$M_R/$file")}
close (D_m);close(D_c);

# ---------------- Merge and summarize -----------------------------
opendir(D_c,"$CAT_R");opendir(D_a,"$A_R/$FILE");opendir(D_m,"$M_R");
	my@file=grep{/^\w/}readdir D_c;
		foreach my$file(@file){open(IN,"$CAT_R/$file");open(ALL,">>$CAT_R/START\_$FILE");select(ALL);
		while(<IN>){print "$file\t$_"}close(IN);close(ALL);
		copy("$CAT_R/START\_$FILE","$A_R/$FILE/START\_$FILE");
		copy("$CAT_R/START\_$FILE","$M_R/$FILE")}
close(D_c);close(D_a);close(D_m);

# ----------------- Summarize and delete ----------------------------
opendir (D_c,"$CAT_R");opendir(D_a,"$A_R/$FILE");

open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
print "\nPossible ATG extensions\nATG extensions and the rest must add up to START\_$FILE total (see below) and m\_$FILE (see above)\n";close(O);

	my@file=grep{/^\w/}readdir D_c;
		foreach my$file(@file){open(IN,"$CAT_R/$file");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
		while(<IN>){}print "$file,$.\n";close(IN);close(Z);unlink("$CAT_R/$file")}
close(D_c);close(D_a);

