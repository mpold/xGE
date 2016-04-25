# All non-SNV adjusted for the left- and right-match
# ONE matching base left on the left side, and none on the right side if Ref and Alt both contain more than 1 base
# if no match found on the left side then nothing done on the left side
# prefixes describing specific variations in LEFT- and RIGHT TRIM are provided in the body of the script

# FINAL OUTPUT: goes into next step: $FILE
# TOTAL of ELEMENTS: 15 (see below)
# [0] - VAR_TYPE:		variation type on genomic DNA
# [1] - CHROM:			chromosome - e.g. 1 (.vcf-s chr is nipped)
# [2] - LEFT_CX_POS:	Corrected start of variation, post LEFT-TRIM = chromosome position + LEFT-MATCH; chromosome position is the corrected value from x_0_M1.pl (upstream script)
# [3] - ID:				period (.) or smth else (comes from .vcf)
# [4] - RL_REF:			RIGHT- and LEFT-trimmed REF
# [5] - RL_ALT:			RIGHT- and LEFT-trimmed ALT
# [6] - REV_RL_ALT:		reversed RIGHT- and LEFT-trimmed (TRIMMED) ALT to determine inversions
# [7] - L_LEN:			LEFT-TRIM length (number of bases trimmed on the left side); this valeue is added to [1], see the very bottom of this script for expanded explanation
# [8] - C_RL_REF:		cut and RIGHT- and LEFT-trimmed (TRIMMED) REF of insdel (see above $CUT), if empty then INSERTION
# [9] - C_RL_ALT:		cut and RIGHT- and LEFT-trimmed (TRIMMED) ALT of insdel (see above $CUT), if empty then DELETION, if both [0] and [1] are empty then variation other than insdel
# [10] - LEN_RL_REF:	length of UNCUT, TRIMMED REF
# [11] - LEN_RL_ALT:	length of UNCUT, TRIMMED ALT
# [12] - CX_POS:		input chromosome position (.vcf-s chr nipped), adjusted for 0- and 1-reference matching in x_0_M1.pl
# [13] - CODE:			input 0- or 1-reference matching code (2,1,21): 2=1-based,1=0-based,21=both
# [14] - DOM:			Nmer domain wherein the variation is found 

# ARCHIVED FILES:
# 1. T_$FILE - contains all entries that underwent trimming etc by this script
# 2. Zummary_$FILE.csv - line counts from x_0_M1.pl and x_SORT_QC.pl
					## OIF = (q + x1) OIF name
					## q = (M_i_nr_q) + (x_i_nr_q) OIF name
					## pre-trim (no 'R_B_L_' files) = post-trim ('R_B_L_ files')
					## pre-trim = post-trim = M_i_nr_q
					## each individual pre-trim file = appropriate post-trim file
					## $FILE - contains data entering mapping to UTR, exons, introns etc.
# 3. $FILE - output from this step; named as the file that enters the pipeline
				
###############################################################################

#!/usr/bin/perl -w

use Data::Dumper;
use File::Copy;
use Parallel::ForkManager;
use POSIX;
use List::Util qw[min max];

my$I_R='';
my$Temp_R='';
my$A_R='';
my$COM_R='';

# -------- STEP1: evaluate LEFT-side-----------------------
opendir(D_i,"$I_R");opendir(D_t,"$Temp_R");
	my@file=grep{/^\w/}readdir D_i;
		foreach my$file(@file){$FILE=$file;
			open(IN,"$I_R/$file");
			# make the rest of the files:
			open(BAD,">>$Temp_R/x2\_$file");close(BAD);
			open(LMNV,">>$Temp_R/LMNV\_$file");close(LMNV);
			open(LSNV,">>$Temp_R/LSNV\_$file");close(LSNV);
			open(mmnv,">>$Temp_R/mmnv\_$file");close(mmnv);
			open(snv,">>$Temp_R/snv\_$file");close(snv);
			open(LINS,">>$Temp_R/LINS\_$file");close(LINS);
			open(LDEL,">>$Temp_R/LDEL\_$file");close(LDEL);
			open(mins,">>$Temp_R/mins\_$file");close(mins);
			open(mdel,">>$Temp_R/mdel\_$file");close(mdel);
			
				while(<IN>){
					s/^chrX/chr23/;			# reverse this when trimming and qc is done!
					s/^chrY/chr24/;			# reverse this when trimming and qc is done!
					# s/^chrM/chr25/;		# FUTURE: add mitochondrial reference sequence to human genome source file
			
					my@VCF=split/chr|\t/;
					my$Ref=$VCF[4];my$Alt=$VCF[5];
					my$Ref_L=length($VCF[4]);my$Alt_L=length($VCF[5]);

# ------------ EXCLUDE POOR QUALITY VARIANTS -------------------------------------
# no difference between Ref and Alt
# chromosomes described as no 1-22, X or y
# chromosomes containing positions larger than the largest possible on Chromosome 1
# no PASS variants
# variants containing bases other than A, C, G, T, U, X, N

my$CHR=$VCF[1];
my$W='\D';
my$L='\W';								# if letters other than X or Y are present in original input
my$NON_NUMB=$_ if/^chr$W/;				# FUTURE: make this condition more strict
my$NON_LETT=$_ if 
($Ref=~/$L/||$Alt=~/$L/);				# no no-letter characters in Ref or Alt
my$COMMENT=$_ if/^\#/;					# .vcf column header and comment lines

my$UN_VCF= 								# describes .vcf format requirements; if not met -> poor quality entry
($VCF[2]<250000000&&					# FUTURE: make the '<' condition chromosome specific 
$VCF[3]=~/(\.|\w)/&& 					# Current version considers only the rough length of chromosome 1 (largest chromosome)
$Ref=~/(A|C|G|T|U|N|X)+/&&
$Alt=~/(A|C|G|T|U|N|X)+/&&
($VCF[6]==1||$VCF[6]==2||$VCF[6]==21)&&
($VCF[7]eq'L'||$VCF[7]eq'M'||$VCF[7]eq'S')&&
$VCF[8]<250000000&&
$VCF[9]=~/(A|C|G|T|U|N|X)+/&&
$CHR<25)&& 								# FUTURE: $CHR<= 25 capable to handle mitochondrial genome as well
!$NON_NUMB&&
!$NON_LETT&&
!$COMMENT;

if($ZERO||$NO_VARIATION||!$UN_VCF){ # exclude entries of poor quality, or no difference between Alt and Ref
	open(BAD,">>$Temp_R/x2\_$file");select(BAD);print;close(BAD)}
		else{ # SORT input into left match/mismatch and length and non-length variation
			my$INS=($Ref_L<$Alt_L)&&($Ref_L ne 0);  # ..in parentheses excludes Ref = '' -> downstream steps don't work if TRUE
			my$DEL=($Ref_L>$Alt_L)&&($Alt_L ne 0);  # ..in parentheses excludes Alt = '' -> downstream steps don't work if TRUE			
			my$MNV=($Ref_L==$Alt_L)&&$Ref_L>1;
			my$SNV=($Ref_L==$Alt_L)&&$Ref_L==1;

			# -------------- LEFT MISMATCH -----------------------------
			my$LEFT=(substr($Ref,0,1) eq substr($Alt,0,1)); # left match
			my$left=(substr($Ref,0,1) ne substr($Alt,0,1)); # left mismatch
		
				if($LEFT&&$MNV){open(LMNV,">>$Temp_R/LMNV\_$file");select(LMNV);print;close(LMNV);		# left- AND right-trim
				}
				elsif($LEFT&&$SNV){open(LSNV,">>$Temp_R/LSNV\_$file");select(LSNV);print;close(LSNV);
				}
				elsif($left&&$MNV){open(mmnv,">>$Temp_R/mmnv\_$file");select(mmnv);print;close(mmnv);	# right-trim
				}
				elsif($left&&$SNV){open(snv,">>$Temp_R/snv\_$file");select(snv);print;close(snv);
				}
				elsif($LEFT&&$INS){open(LINS,">>$Temp_R/LINS\_$file");select(LINS);print;close(LINS);	# left- AND right-trim
				}
				elsif($LEFT&&$DEL){open(LDEL,">>$Temp_R/LDEL\_$file");select(LDEL);print;close(LDEL);	# left- AND right-trim
				}
				elsif($left&&$INS){open(mins,">>$Temp_R/mins\_$file");select(mins);print;close(mins);	# right-trim
				}
				elsif($left&&$DEL){open(mdel,">>$Temp_R/mdel\_$file");select(mdel);print;close(mdel);	# right-trim
				}}}
	close(IN); # Close the original input file
}close(D_i);close(D_t);


# ------------- Summarize folder -------------------
opendir(D_a,"$S_R");opendir(D_t,"$Temp_R");

open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
print"\nLEFT_TRIM variations: L-variations except LSNV\_$FILE\n";
# print "RIGHT_TRIM variations: all variation except LSNV\_ and snv\_$FILE\n";
close(O);

	my@f=grep{/^\w/}readdir D_t;
		foreach my$f(@f){open(I,"$Temp_R/$f");open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
		while(<I>){}print"$f,$.\n";close(I);close(O)}
close(D_t);close(D_a);

# STEP2: LEFT-trim - greater than 1 base LEFT-match is reduced to a single base LEFT-match in length variants, and to 0 match in MNV
# Left-trimming QC-ed: variants with trimmed left side and adjusted coordinates produced the same results as unadjusted and -trimmed (01-18-2016)
opendir(D_t,"$Temp_R");
	my@file=grep{/^(LMNV\_|LINS\_|LDEL\_)/}readdir D_t;	# see above for the prefixes for LEFT-trimming
		foreach my$file(@file){open(IN,"$Temp_R/$file");
		
# make output files:
open(O,">>$Temp_R/L\_$file\.csv");close(O);
		
		while(<IN>){
			my@V=split/\t/;
			my$Ref=$V[3];my$Alt=$V[4];
			my$Ref_L=length($V[3]);my$Alt_L=length($V[4]);

			my$INS=($Ref_L<$Alt_L)&&($Ref_L ne 0);	# the part in parentheses excludes entries with Ref = '' -> downstream steps don't work
			my$DEL=($Ref_L>$Alt_L)&&($Alt_L ne 0);	# the part in parentheses excludes entries with Alt = '' -> downstream steps don't work
			my$MNV=$Ref_L==$Alt_L;					# MNV with Ref and Alt first base matching
			# ----------------------------------------------------------------------
				if($INS||$DEL){open(INSDEL,">>$Temp_R/L\_$file\.csv");select(INSDEL);
					my$TOP=max($Ref_L,$Alt_L);

					# --------- LEFT_MATCH and _TRIM --------------------------
					foreach(my$L_M=1;$L_M<=$TOP;$L_M++){
						my$LEFT_MATCH=(
						(substr($V[3],0,$L_M) eq substr($V[4],0,$L_M))&&
						(substr($V[3],0,($L_M+1)) ne substr($V[4],0,($L_M+1))));

						if ($LEFT_MATCH){
							my$L_Trim=$L_M-1;
							my$Cx_GRCh=$V[1]+$L_Trim;	# $V[1] - position [potentially] corrected by x_0_M1.pl
							my$Cx_Ref=substr($Ref,$L_Trim,($Ref_L-$L_Trim));
							my$Cx_Alt=substr($Alt,$L_Trim,($Alt_L-$L_Trim));

print"$L_Trim,$L_M,$Cx_Ref,$Cx_Alt,$Cx_GRCh,\($V[3]\),\($V[4]\),$V[0],$V[1],\($V[2]\),$V[5],\($V[6]\),$V[7],\($V[8]\)\n";
							}}close(INSDEL);
							}else{
							
							open(MNV,">>$Temp_R/L\_$file\.csv");select(MNV); # else can be used here because nothing else is left to process
							my$TOP=max($Ref_L,$Alt_L);

							# -------- LEFT_MATCH and _TRIM -----------------------
							foreach(my$L_M=1;$L_M<=$TOP;$L_M++){
							my$LEFT_MATCH=(
							(substr($V[3],0,$L_M) eq substr($V[4],0,$L_M))&&
							(substr($V[3],0,($L_M+1)) ne substr($V[4],0,($L_M+1))));

							if($LEFT_MATCH){
								my$L_Trim=$L_M; # $L_M is greater by one for MNV because Zero-left match is needed 
												# most MNV turn into SNV
								my$Cx_GRCh=$V[1]+$L_Trim;	# corrected chromosome position
								my$Cx_Ref=substr($Ref,$L_Trim,($Ref_L - $L_Trim));
								my$Cx_Alt=substr($Alt,$L_Trim,($Alt_L - $L_Trim));
print"$L_Trim,$L_M,$Cx_Ref,$Cx_Alt,$Cx_GRCh,\($V[3]\),\($V[4]\),$V[0],$V[1],\($V[2]\),$V[5],\($V[6]\),$V[7],\($V[8]\)\n";
						}}close(MNV)}}close(IN)}
close (D_t);
# -----------------------------------------------------------------------------------------------------------------------
# STEP3: format -- equalize the number of elements -- the non-left trim files the same as the left-trim files (see above)
# Otherwise, the downstream steps become difficult

# 		LMNV # left- AND right-trim
# LSNV
# mmnv # right-trim
# snv
# 		LINS # left- AND right-trim
# 		LDEL # left- AND right-trim
# mins # right-trim
# mdel # right-trim

opendir(D_t,"$Temp_R");
	my@file=grep{/^(LSNV\_|mmnv\_|snv\_|mins\_|mdel\_)/}readdir D_t;	# see above for the prefixes for LEFT-trimming
	foreach my$file(@file){open(IN,"$Temp_R/$file");open(OUT,">>$Temp_R/$file\.csv");select(OUT);
	while(<IN>){my@V=split/\t/;
	print"0,0,$V[3],$V[4],$V[1],\($V[3]\),\($V[4]\),$V[0],$V[1],\($V[2]\),$V[5],\($V[6]\),$V[7],\($V[8]\)\n";
	}close(IN);close(OUT)}
close(D_t);

# ----------- STEP4: Block first base of Ref and Alt to right-trim -------------------------------
# Delete all B_ files at the end
opendir(D_t,"$Temp_R");
	my@file=grep{/\.csv/}readdir D_t;	# see above for the prefixes for LEFT-trimming
	foreach my$file(@file){open(IN,"$Temp_R/$file");open(BLOCK,">>$Temp_R/B\_$file");select(BLOCK);
		while(<IN>){
		# Block first base of Ref
		s/,A/,I/; s/,C/,J/; s/,G/,K/; s/,T/,L/; s/,U/,M/; s/,N/,N/; s/,X/,O/;
		# Block first base of Alt
		s/,A/,I/; s/,C/,J/; s/,G/,K/; s/,T/,L/; s/,U/,M/; s/,N/,N/; s/,X/,O/;
		print}close(IN);close(BLOCK);unlink("$Temp_R/B\_Zummary.csv")}
close(D_t);


# ----------- STEP5: Right-trim wherever necessary -----------------------------------
# THE Right-trim files: 6 IN TOTAL
# B_L_LMNV # left- AND right-trim
# B_mmnv # right-trim
# B_L_LINS # left- AND right-trim
# B_L_LDEL # left- AND right-trim
# B_mins # right-trim
# B_mdel # right-trim

opendir(D_t, "$Temp_R");
	my@file=grep{/^(B\_)/}readdir D_t;	# see above for the prefixes for LEFT-trimming
		foreach my $file (@file){open(IN,"$Temp_R/$file");
		# make the rest of the files:
		open(O,">>$Temp_R/R\_$file");close(O);
			while(<IN>){
				my@V_Cx=split/\,/;
				my$Ref_Cx=$V_Cx[2];my$Alt_Cx=$V_Cx[3];
				my$Ref_L_Cx=length($V_Cx[2]);my$Alt_L_Cx=length($V_Cx[3]);
				my$Ref_z=substr($V_Cx[2],-1,1); # last Ref base
				my$Alt_z=substr($V_Cx[3],-1,1); # last Alt base
########## DEFINE VARIATIONS: redundant but good for qc, never know what these files contain
my$INS=$Ref_L_Cx<$Alt_L_Cx;
my$DEL=$Ref_L_Cx>$Alt_L_Cx;
my$MNV_SNV=$Ref_L_Cx==$Alt_L_Cx;
my$LIM=(min($Ref_L_Cx,$Alt_L_Cx))-1;
				if($INS||$DEL||$MNV_SNV){open(VAR,">>$Temp_R/R\_$file");select(VAR);	# opens all appropriate files
					# NO RIGHT-TRIM needed!
					if($Ref_z ne $Alt_z){print"0,$Ref_L_Cx,$Alt_L_Cx,$Ref_Cx,$Alt_Cx,$_"}
					else # takes on files in which right-match and -trim is needed (Ref_z eq Alt_z)
						{
						foreach(my$i=-1;$i>=-$LIM;$i--){
							if(((substr($V_Cx[2],$i,-$i) eq substr($V_Cx[3],$i,-$i))&&
								(substr($V_Cx[2],($i-1),-($i-1)) ne substr($V_Cx[3],($i-1),-($i-1))))){
									my$Ref_mism=substr($V_Cx[2],0,($Ref_L_Cx + $i));
									my$Alt_mism=substr($V_Cx[3],0,($Alt_L_Cx + $i));
									print"$i,$Ref_L_Cx,$Alt_L_Cx,$Ref_mism,$Alt_mism,$_"}}}}
					close(VAR)}close(IN)}
close(D_t);

# ------ STEP7: Unlink pre-cursor files ---------------
opendir(D_x,"$Temp_R");
	my@file=grep{!/^(R\_B\_|Zummary|x1\_)/}readdir D_x;
		foreach my$file(@file){unlink("$Temp_R/$file")}
close(D_x);

# ------------- STEP8: Re-format ------------------------------
opendir(D_i,"$I_R");opendir(D_t,"$Temp_R");opendir (D_a,"$A_R");
	my@file=grep{/^R\_/}readdir D_t;
		foreach my$file(@file){open(IN,"$Temp_R/$file");open(RE,">>$Temp_R/T\_$FILE\.csv");select(RE);
			while(<IN>){
			s/,I/,A/g;s/,J/,C/g;s/,K/,G/g;s/,L/,T/g;s/,M/,U/;s/,N/,N/;s/,O/,X/;
			s/,\(/,/g;s/\),/,/g;s/\)$//;
			s/chr23/chrX/;s/chr24/chrY/; # restore X and Y

			my@Fin_TRIM=split/\,/;
			my$Fin_R_L=length($Fin_TRIM[3]);my$Fin_A_L=length($Fin_TRIM[4]);
			my$CUT_R=substr($Fin_TRIM[3],1,($Fin_R_L-1));
			my$CUT_A=substr($Fin_TRIM[4],1,($Fin_A_L-1));

			my$CUT= # report format (cut) deletions and insertions (either Ref or Alt is one base, post trimming (see above)
			(($Fin_R_L ne $Fin_A_L)&&
			(substr($Fin_TRIM[3],0,1) eq substr($Fin_TRIM[4],0,1)));

			if($CUT){print"$CUT_R,$CUT_A,$Fin_R_L,$Fin_A_L,$_"}
			else{print",,$Fin_R_L,$Fin_A_L,$_"; # variants other than insertions and deletions
			}}close(IN);close(RE)}

# ---------- FINALIZE OUTPUT (see the top of this file for details) ---------------------------
# Variation type as part of final output is added in a downstream step
open(IN,"$Temp_R/T\_$FILE\.csv");open(RE,">>$I_R/F\_$FILE");select(RE);
	while(<IN>){my@FIN=split/\,/;
		print"$FIN[16]\t$FIN[13]\t$FIN[18]\t$FIN[7]\t$FIN[8]\t";
		print scalar reverse"$FIN[8]";print"\t";	# reversed ALT, in order to determine inversions
		print"$FIN[9]\t$FIN[0]\t$FIN[1]\t$FIN[2]\t$FIN[3]\t$FIN[17]\t$FIN[19]\t$FIN[22]"}
	close(IN);close(RE);copy("$Temp_R/T\_$FILE\.csv","$A_R/$FILE/T\_$FILE\.csv");
close(D_t);
	unlink ("$I_R/$FILE"); # unlink the input of this script, it already exist in the archive folder
close(D_i);close(D_a);


# ----------- STEP9: Summarize folder ----------
opendir(D_a,"$A_R");opendir(D_t,"$Temp_R");
	unlink("$Temp_R/R\_B\_Zummary\_$FILE\.csv");

open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
print"\nRIGHT_TRIM variations: all input variation except LSNV\_ and snv\_$FILE \(see above\)\n";
print"SUM of R\_B_files must be total of $FILE \(see below\)\n";
close(OUT);

	my@file=grep{/^R\_/}readdir D_t;
		foreach my$file(@file){open(IN,"$Temp_R/$file");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
		while(<IN>){}print"$file,$.\n";close(IN);close(OUT)}

# --------- STEP10: clean up the folders --------------------------------------------
	system("$COM_R/CLEAN_UP.pl");	 # alter this as needed in the production version
close(D_t);close(D_a);


# --------- STEP10: determine variation types on genomic DNA level -------------------
opendir(D_i,"$I_R");opendir(D_t,"$Temp_R");
	my@file=grep{/^F\_/}readdir D_i;
		foreach my$file(@file){

		open(IN,"$I_R/$file");
		# ----- make the rest of the files ---------
		open(S,">>$Temp_R/SNV\_$file\.csv");close(S);
		open(V,">>$Temp_R/INV\_$file\.csv");close(V);
		open(M,">>$Temp_R/MNV\_$file\.csv");close(M);
		open(U,">>$Temp_R/DUP\_$file\.csv");close(U);
		open(I,">>$Temp_R/INS\_$file\.csv");close(I);
		open(L,">>$Temp_R/DEL\_$file\.csv");close(L);
		open(CLX,">>$Temp_R/COMPLEX\_$file\.csv");close(CLX);
		
		while(<IN>){s/^chr//g;s/\t/,/g;

		my@F =split/\,/;

		my$REF=$F[3];my$REF_L=length($F[3]);my$REF_1=substr($F[3],0,1);
		my$ALT=$F[4];my$ALT_L=length($F[4]);my$ALT_1=substr($F[4],0,1);
		my$ALT_INV=$F[5];	# reversed ALT to determine inversions
		my$REF_T=$F[7];my$REF_T_L=length($REF_T);
		my$ALT_T=$F[8];my$ALT_T_L=length($ALT_T);
		my$NMER=$F[13];
		my$DIF=$F[1]-$F[11]; # $DIF = J1 in Excel formula below

		my$sub_SNV=($REF_L==$ALT_L)&&$REF_L==1;
		my$INV=$REF eq $ALT_INV;
		my$sub_MNV=($REF_L==$ALT_L)&&$REF_L>1&&!$INV;

		# EXCEL CONTROL FORMULA: =IF(MID(X1,(100+J1),LEN(L1))=MID(L1,1,LEN(L1)),TRUE,FALSE)
		# 0-and 1-based condition removed in 02162016 because the 1-shift for 0-based carried in x_adj_NMER.pl
		my$DUP=
		(($ALT_T eq (substr($NMER,(100+$DIF-$ALT_T_L),$ALT_T_L)))&&
		($REF_L<$ALT_L&&$ALT_T_L>0));

		my$INS=($REF_L<$ALT_L)&&!$DUP&&$REF_L==1&&$REF_1 eq $ALT_1;
		my$DEL=($REF_L>$ALT_L)&&$ALT_L==1&&$REF_1 eq $ALT_1;

		if ($sub_SNV){open(S,">>$Temp_R/SNV\_$file\.csv");select(S);print"SNV,$_";close(S);
			open(O,">>$I_R/$FILE");select(O);print"SNV,$_";close(O)}
		elsif($INV){open(V,">>$Temp_R/INV\_$file\.csv");select(V);print"inv,$_";close(V);
			open(O,">>$I_R/$FILE");select(O);print"inv,$_";close(O)}
		elsif($sub_MNV){open(M,">>$Temp_R/MNV\_$file\.csv");select(M);print"MNV,$_";close(M);
			open(O,">>$I_R/$FILE");select(O);print"MNV,$_";close(O)}
		elsif($DUP){open(U,">>$Temp_R/DUP\_$file\.csv");select(U);print"dup,$_";close(U);
			open(O,">>$I_R/$FILE");select(O);print"dup,$_";close(O)}	
		elsif($INS){open(I,">>$Temp_R/INS\_$file\.csv");select(I);print"ins,$_";close(I);
			open(O,">>$I_R/$FILE");select(O);print"ins,$_";close(O)}
		elsif($DEL){open(L,">>$Temp_R/DEL\_$file\.csv");select(L);print"del,$_";close(L);
			open(O,">>$I_R/$FILE");select(O);print"del,$_";close(O)}
		else{open(CLX,">>$Temp_R/COMPLEX\_$file\.csv");select(CLX);print"complex,$_";close(CLX);	# in case the variations as defined above produce no results
		open(O,">>$I_R/$FILE");select(O);print"complex,$_";close(O)}}
	close(IN)}
close(D_i);close(D_t);

# ------------- STEP11: Summarize folder ---------------------
opendir(D_a,"$A_R");opendir(D_t,"$Temp_R");opendir(D_i,"$I_R");
	
open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
print"\n$FILE and F\_$FILE must be equal, and equal to M\_i\_nr\_q\_$FILE \(see the top of this file\)\n";
print"M\_i\_nr\_q\_$FILE contains all upstream QC\-pass entries\n";
print"Variation types below must add up to $FILE\n";
close(OUT);
	
	my@file=grep{/^\w/}readdir D_t;
		foreach my$file(@file){open(IN,"$Temp_R/$file");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
			while(<IN>){}print"$file,$.\n";close(IN);close(OUT)}
	my@f=grep{/^\w/}readdir D_i;
		foreach my $f(@f){open(IN,"$I_R/$f");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
			while (<IN>){}print "$f,$.\n";close(IN);close(OUT)}
##############################
# STEP11: clean up the folders
	system("$COM_R/CLEAN_UP.pl");	 # alter this as needed in the production version
	copy("$I_R/$FILE","$A_R/$FILE");
	unlink("$I_R/F\_$FILE");
close(D_t);close(D_a);close(D_i);