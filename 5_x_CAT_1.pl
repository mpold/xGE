# x_CAT_1.pl determines possible UTR5 and UTR3 variants affecting ATG and STOP, respectively
# It also determines possible coding region variants affecting ATG and STOP
# Also, distances from cdsSTART and cdsEND are determined in the coding region files
# Also, this the first step to determine the entries that go into downstream TRANSLATION step:
		# 'pC' and 'mC' prefix files: plus and minus coding region
		# 'atg_u5pC and 'stop_u3pC': utr5 PLUS variants that may affect ATG, and utr3 MINUS variants that may affect STOP
# 
# ARCHIVE: cat_1_file - contains the merged results of this step; if needed can be split into sub-categories as during the execution of this script_name
#          # in addition to variables described below, the file name of each split is added in front of every entry in the merger

# Three (4) elements are added to each entry in the current step: 
# 1 tag (file prefix (the leftmost code, and 3 computed values (three rightmost values))
# FILE PREFIXES:...........................................................VARIABLES ADDED TO INPUT
# atg_u5pC  - utr5 variants affecting ATG, PLUS............................U5_P,$ATG_P,,
# u5p       - all other utr5, PLUS.........................................U5_P,$ATG_P,,
# u3p       - all utr3, PLUS...............................................U3_P,,$STOP_P,
# 
# stop_u3mC - utr3 variants affecting STOP, MINUS..........................U3_M,$STOP_M,,
# u3m       - all other utr3, MINUS........................................U3_M,$STOP_M,,
# u5m       - all utr5, MINUS..............................................U5_M,,$ATG_M,
# 
# ncRNA_P   - PLUS non-coding rna..........................................ncRNA_P,,,
# ncRNA_M   - MINUS non-coding rna.........................................ncRNA_M,,,
# 
# atg_pC    - ATG affectors, PLUS coding region............................CDS_P,$ATG_P,$STOP_P,
# stop_pC   - STOP affectrs, PLUS coding region............................CDS_P,$ATG_P,$STOP_P,
# pC        - PLUS coding region, other than ATG and STOP affectors........CDS_P,$ATG_P,$STOP_P,
# 
# atg_mC    - ATG affectors, MINUS coding region...........................CDS_M,$STOP_M,$ATG_M,
# stop_mC   - STOP affectors, MINUS coding region..........................CDS_M,$STOP_M,$ATG_M,
# mC        - MINUS coding region, other than ATG and STOP affectors.......CDS_M,$STOP_M,$ATG_M, 
# 
# int       - entries that did not fit into any of the above categories....intergenic,,,
#             STOP, if you see anything here, upstream steps may not work properly

# NB! - taking account the 0-based system is not necessary because only the UCSC GB cdsSTART and cdsEND matter

# OUTPUT ELEMENTS: 19 + UCSC exome descriptions
# [0] - TAG_1: transcript mapping file prefix
# [1] - MAP: computed transcript mapping values as shown above (on the right side)
# [2] - AP_SM: distance from either plus strand ATG or minus strand STOP
# [3] - AM_SP: distance from either minus strand ATG or plus strand STOP
# [4] - [18] - the [0] - [14] from x_SORT.pl; all elements now tab-separated followed by..
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
		foreach my$file(@file){$FILE = $file;
# STEP1: map all variants from that qualify from the previous steps: 
			open(IN,"$M_R/$file");
				while(<IN>){
no warnings;			
				# DO X->23 and Y->24 formatting, needed for downstream 1-to-24-loops
				s/\tX\t/\t23\t/;s/\tY\t/\t24\t/;
				s/\tchrX\t/\t23\t/;s/\tchrY\t/\t24\t/;
				s/\tchr/\t/;

				my @U = split /\t|,/; # first comma denotes the start of UCSC GB data

# PART 1: map out the NM and NR EXON/INTRON placements of variants
# add 11 to all UCSC GB elements
# U[2] - variant
# U[19] - TxStart
# U[20] - TxEnd
# U[21] - cdsStart
# U[22] - cdsEnd

my$Tx_L=$U[20]-$U[19];
my$CDS_L=$U[22]-$U[21];
my$U_5P=(($U[2]>$U[19])&&($U[2]<$U[21])&&$U[18]=~/\+/)if$U[16]=~/NM\_/;
my$U_3P=(($U[2]>$U[22])&&($U[2]<($U[20]+1))&&$U[18]=~/\+/)if$U[16]=~/NM\_/;
my$U_5M=(($U[2]>$U[22])&&($U[2]<($U[20]+1))&&$U[18]=~/\-/)if$U[16]=~/NM\_/;
my$U_3M=(($U[2]>$U[19])&&($U[2]<($U[21]))&&$U[18]=~/\-/)if$U[16]=~/NM\_/;
my$ncRNA_P=$U[18]=~/\+/ if$U[16]=~/NR\_/;
my$ncRNA_M=$U[18]=~/\-/ if$U[16]=~/NR\_/;
my$cds_P=(($U[2]>=$U[21]&&$U[2]<=$U[22])&&$U[18]=~/\+/)if$U[16]=~/NM\_/;
my$cds_M=(($U[2]>=$U[21]&&$U[2]<=$U[22])&&$U[18]=~/\-/)if$U[16]=~/NM\_/;

# START CODON:
# ATG_P = 0 means that ATG positions are 1,2,3 (ATG/CDS start on PLUS is at 1-position)
# ATG_P = -1 means that the A of ATG is the SECOND position from VARIANT position
# ATG_P should not be negative on CDS variants, only UTR5_PLUS produces negative vales
my$ATG_P=($U[2]-$U[21])if$U[18]=~/\+/;		# PLUS strand: if negative then upstream (TO); if positive then downstream (FROM)

# ATG_M = -2 means that CAT positions are -2,-1,0 (ATG/CDS start on MINUS is at 0-position)
# ATG_M should not be positive on CDS variants, only UTR5_MINUS produces positive values
my$ATG_M=($U[2]-$U[22])if$U[18]=~/\-/;		# MINUS strand: if negative then downstream (FROM); if positive the upstream (TO) 
													# They matter in UTR5 and CDS
													# Don't matter in ncRNA, don't even write them out
# STOP CODON:
# STOP_P = -2 means that FIRST position of STOP is variant
# STOP_P = 0 means that LAST position of STOP is variant -> CDS end is 0-based
my$STOP_P=($U[2]-$U[22])if$U[18]=~/\+/;	# PLUS: if negative then upstream (TO); if positive then downstream (FROM)

# STOP_M = 0 means the first downstream base from STOP -> CDS end is 1-based
my$STOP_M=($U[2]-$U[21])if$U[18]=~/\-/;	# MINUS: if negative then downstream (FROM); if positive the upstream (TO)
													# They matter in UTR3 and CDS
													# They don't matter in ncRNA, don't write them out
my$Ref_L=length($U[4]);my$Alt_L=length($U[5]);

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
# 
# 
# FILE PREFIXES:...........................................................VARIABLES ADDED TO INPUT
# atg_u5pC  - utr5 variants affecting ATG, PLUS............................U5_P,$ATG_P,,
# u5p       - all other utr5, PLUS.........................................U5_P,$ATG_P,,
# u3p       - all utr3, PLUS...............................................U3_P,,$STOP_P,
#
# stop_u3mC - utr3 variants affecting STOP, MINUS..........................U3_M,$STOP_M,,
# u3m       - all other utr3, MINUS........................................U3_M,$STOP_M,,
# u5m       - all utr5, MINUS..............................................U5_M,,$ATG_M,
#
# ncRNA_P   - PLUS non-coding rna..........................................ncRNA_P,,,
# ncRNA_M   - MINUS non-coding rna.........................................ncRNA_M,,,
#
# atg_pC    - ATG affectors, PLUS coding region............................CDS_P,$ATG_P,$STOP_P,
# stop_pC   - STOP affectrs, PLUS coding region............................CDS_P,$ATG_P,$STOP_P,
# pC        - PLUS coding region, other than ATG and STOP affectors........CDS_P,$ATG_P,$STOP_P,
#
# atg_mC    - ATG affectors, MINUS coding region...........................CDS_M,$STOP_M,$ATG_M,
# stop_mC   - STOP affectrs, MINUS coding region...........................CDS_M,$STOP_M,$ATG_M,
# mC        - MINUS coding region, other than ATG and STOP affectors.......CDS_M,$STOP_M,$ATG_M, 
#
# int       - entries that did not fit into any of the above categories....intergenic,,,
#             STOP, if you see anything here, upstream steps may not work properly

					
# -------------- Make all output files: file prefixes make up the array -----------------------------------------------------------
my@OUT_F=('atg_u5pC','u5p','u3p','stop_u3mC','u3m','u5m','ncRNA_P','ncRNA_M','atg_pC','stop_pC','pC','atg_mC','stop_mC','mC','int');
foreach my$OUT_F(@OUT_F){open(O_F,">>$CAT_R/$OUT_F\_$file");close(O_F)}

				# DETERMINE variants originating from UTR and affecting CDS: DELETIONS, MNV, INV, and COMPLEX
				# They will go into TRANSLATION along with CDS variants

				# UTR5 PLUS
				if($U_5P){																		# deletions close to ATG may affect coding sequence
					if(($Ref_L-1)>abs($ATG_P)){open(O,">>$CAT_R/atg\_u5pC\_$file");select(O);	# values should range from NEGATIVE to ZERO
					print"U5\_P\t$ATG_P\t\t$_";close(O)}
						else{open(O,">>$CAT_R/u5p\_$file");select(O);print "U5\_P\t$ATG_P\t\t$_";close(O)}}
				# UTR5 MINUS
				elsif($U_5M){open(O,">>$CAT_R/u5m\_$file");select(O);print "U5\_M\t\t$ATG_M\t$_";close(O)}
				# UTR3 PLUS
				elsif($U_3P){open(O,">>$CAT_R/u3p\_$file");select(O);print "U3\_P\t\t$STOP_P\t$_";close(O)}
				# UTR3 MINUS
				elsif($U_3M){																	# deletions close to STOP may affect coding sequence
					if (($Ref_L-1)>abs($STOP_M)){open(O,">>$CAT_R/stop_u3mC\_$file");select(O);	# values should range from NEGATIVE to ZERO
					print "U3\_M\t$STOP_M\t\t$_";close(O)}
						else{open(O,">>$CAT_R/u3m\_$file");select(O);print "U3\_M\t$STOP_M\t\t$_";close(O)}}
				# non-coding RNA PLUS
				elsif($ncRNA_P){open(O,">>$CAT_R/ncRNA\_P\_$file");select(O);print "R_P\t\t\t$_";close(O)}
				# non-coding RNA MINUS
				elsif($ncRNA_M){open(O,">>$CAT_R/ncRNA\_M\_$file");select(O);print "R_M\t\t\t$_";close(O)}
				# PLUS coding region including intron variants; the right side of || computes SNV
				elsif($cds_P){
					if((($Ref_L)>=abs($ATG_P))||(($ATG_P<=3)&&($ATG_P>=1))){open(O,">>$CAT_R/atg_pC\_$file");select(O);
					print "CDS\_P\t$ATG_P\t$STOP_P\t$_";close(O)}
						elsif(($Ref_L)>=abs($STOP_P)||(($STOP_P<=0)&&($STOP_P>=-2))){open(O,">>$CAT_R/stop_pC\_$file");select(O);
						print "CDS\_P\t$ATG_P\t$STOP_P\t$_";close(O)}
							else{open(O,">>$CAT_R/pC\_$file");select(O);print "CDS\_P\t$ATG_P\t$STOP_P\t$_";close(O)}}
				# MINUS coding region including intron variants; the right side of || computes SNV
				elsif($cds_M){
					if((($Ref_L)>=abs($ATG_M))||(($ATG_M>=-2)&&($ATG_M<=0))){open(O,">>$CAT_R/atg_mC\_$file");select(O);
					print"CDS\_M\t$STOP_M\t$ATG_M\t$_";close(O)}
						elsif(($Ref_L)>=abs($STOP_M)||(($STOP_M<=3)&&($STOP_M>=1))){open(O,">>$CAT_R/stop_mC\_$file");select(O);
						print"CDS\_M\t$STOP_M\t$ATG_M\t$_";close(O)}
							else{open(O,">>$CAT_R/mC\_$file");select(O);print "CDS\_M\t$STOP_M\t$ATG_M\t$_";close(O)}}
				# The rest: if all is well then nothing should fall into this category
				else{open(O,">>$CAT_R/int\_$file");select(O);print"intergenic\t";close(O)} # theoretically at this point should not pop up at all because of x_INPUT.pl
		}close(IN);close(OUT);unlink("$M_R/$file")}
close (D_m);close(D_c);

# ------------------ Merge and summarize --------------------------
opendir(D_c,"$CAT_R");opendir(D_a,"$A_R/$FILE");opendir(D_m,"$M_R");
	my@file=grep{/^\w/}readdir D_c;
		foreach my$file(@file){open(IN,"$CAT_R/$file");open(ALL,">>$CAT_R/cat\_1\_$FILE");select(ALL);
		while(<IN>){print "$file\t$_";}close(IN);close(ALL);
		copy("$CAT_R/cat\_1\_$FILE","$A_R/$FILE/cat\_1\_$FILE");
		copy("$CAT_R/cat\_1\_$FILE","$M_R/$FILE");
	}unlink("$CAT_R/cat\_1\_$FILE");
close(D_c);close(D_a);close(D_m);

# ----------------- Summarize and delete ---------------------------
opendir (D_c,"$CAT_R");opendir(D_a,"$A_R/$FILE");

open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
print "\nUTR5\_CDS\_UTR3 and ncRNA mapping results\n";
print "All mapping categories must add up to m\_$FILE total (see above)\n";close(O);

	my@file=grep{/^\w/}readdir D_c;
	foreach my $file(@file){open(IN,"$CAT_R/$file");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
	while(<IN>){}print "$file,$.\n";close(IN);close(Z);unlink("$CAT_R/$file")}
close(D_c);close(D_a);
