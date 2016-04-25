######################## PLUS STRAND ##################################
# EXON:   +value: ACCEPTOR side: distance from PREVIOUS intron end (Iz)
# EXON:   -value: DONOR side: distance to NEXT intron start (I1)
# INTRON: +value: DONOR side: distance from PREVIOUS exon end (Ez)
# INTRON  -value: ACCEPTOR side: distance to NEXT exon start (E1)
# 
#  PLUS
#  ------->
#  - - - intron - - - ->|- - - - - - EXON - - - - - ->|- - - - intron - - - ->|- - - - - EXON - - - -
#  . . . . . . . . . . .|1 2 3 . vE. . . . . . .-3-2-1|1 2 3 . . vI . . -1-2-3| . . . . . . . . . . .
#  . . . . . . . . . . .|$From . . . . . . . . . . $To|$From . . . . . . . $To| . . . . . . . . . . .
# 
######################## MINUS STRAND #################################
# EXON:   -value: ACCEPTOR side: distance from PREVIOUS intron end (Iz)
# EXON:   +value: DONOR side: distance to NEXT intron start (I1)
# INTRON: -value: DONOR side: distance from PREVIOUS exon end (Ez)
# ITRON:  +value: distance to NEXT exon start (E1)
# 																          MINUS
#                                                                       <-------
#  - - - intron - - - - |<- - - - - - EXON - - - - - -|<- - - - intron - - - -|<- - - - - EXON - -
#  . . . . . . . . . . .|1 2 3 . vE. . . . . . .-3-2-1|1 2 3 . . vI . . -1-2-3| . . . . . . . . . . 
#  . . . . . . . . . . .|$To . . . . . . . . . . $From|$To . . . . . . . $From| . . . . . . . . . .					
# 
#
# OUTPUT FILES: 
# Prefixes: A - all, b-splice junction and conserved intron bases not affected, JIE - exon/intron junction (splice site crossed)
# 			p - PLUS, m - MINUS, e - exon, i - intron, x - possible affected, n - not affected or reconstituted
#
# 1. A\_$FILE:			see the body of this script and Zummary-file for more detail 
# 2. b\_$FILE:				"
# 3. I\_1z\_$FILE:			"
# 4. I\_2y\_$FILE:			"
# 5. JIE\_$FILE:		All JIE entries will go into translation
# 6. JIE_pix\_$FILE:		"
# 7. JIE_pex\_$FILE:		"
# 8. JIE_pni\_$FILE:		"
# 9. JIE_pne\_$FILE:		"
# 10. JIE_mix\_$FILE:		"
# 11. JIE_mex\_$FILE:		"
# 12. JIE_mni\_$FILE:		"
# 13. JIE_mne\_$FILE:		"
# 
# ARCHIVED FILES: SPLICE_$FILE - merger of the above three files
# 
# OUTPUT ELEMENTS: 31 + UCSC exome descriptions
# Elements from this script:
# [0] - SPL: one of the prefixes, except the 'A' from the above output files; values can be both positive and negative (see the top of this script)
# [1] - EI_LEN: $INEX_L, exon or intron length wherein variant is found
# [2] - Ref_EI_2: $AG_O_P or $AG_O_M: reference 2 first/last intron bases (P - plus, M - minus), value present only in JIE-prefix-files
# [3] - Alt_EI_2: $AG_V_P of $AG_V_M: variant 2 first/last splice intron bases (P - plus, M - minus), value present only in JIE-prefix-files
# [4] - J_Ref_7: $J_O_P or $J_O_M: reference splice junction, 7 bases long, conserved intron bases are 2 and 3 (in perl 0-system), value present only in JIE-prefix-files
# [5] - J_Alt_7: $J_V_P or $J_V_M: variant splice junction, 7 bases long, conserved intron bases are 2 and 3 (in perl 0-system), value present only in JIE-prefix-files
# [6] - FROM: see above for explanation to $From and $To; see the body of this script for more detail on $From and $To
# [7] - TO: 
# [8] - EI_id: exon -> E; intron -> I
# [9] - V_EI: the variant exon/intron number
# Elements from the previous script:
# [10] - [30] and all elements of UCSC exome entries

#!/usr/bin/perl -w
use Data::Dumper;
use File::Copy;
use POSIX;
use List::Util qw[min max];
use File::Copy;


my$CAT_R='';
my$M_R='';
my$I_R='';
my$A_R='';

opendir(D_m,"$M_R");opendir (D_cat,"$CAT_R");
	my@file=grep {/^(\w|\d)/}readdir D_m;
		foreach my$file(@file){$FILE = $file;	
# STEP1: determine the exon/intron distribution of variants
			open(IN,"$M_R/$file");open(O,">>$CAT_R/A\_$file");select(O);			
# make additional output files:
my @OUT_F = ('JIE','I\_1z','I\_2y','b');
foreach my $OUT_F(@OUT_F){open(OUT_F,">>$CAT_R/$OUT_F\_$FILE");close(OUT_F)}
					while(<IN>){
# no warnings;
						my @U=split/\,|\t/;
# my $Tx_L = abs($U[25] - $U[26]);
# my $CDS_L = abs($U[27] - $U[28]);
# my $Viant = $U[8];		# position on chromosome
# my $Str = $U[24];
#
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

						# The plus strand constant: computed based on the previous versions of this script
						my$CONSTANT_P=31;	# must by greater by 2 as compared to $Exon_Count (see below)

						# Exon count:
						my$Exon_Count=$U[29];

							foreach(my$i=($CONSTANT_P-1);$i<=(($CONSTANT_P-1)+($Exon_Count-1));$i++){
							my$j=$i+($Exon_Count+1);my$Exon_P=$i-29;my$Exon_M=$Exon_Count-($i-30);

# ----------- PLUS STRAND: the P in the below variables stands for PLUS strand ----------------------
# $i           $j         $i+1          $j+1    
# E1E2->...EyEz->I1I2...IyIz->E1E2...EyEz->I1I2

# ----> PLUS DIRECTION
# ------ Exon ------|______ INTRON______|------ Exon ------
#                 Donor               Acceptor

my$Az_E_P=$U[8]-$U[$i];			# +vale: ACCEPTOR side: distance from PREVIOUS intron end (Iz)
my$D1_E_P=$U[8]-($U[$j]+1);		# -value: DONOR side: distance to NEXT intron start (I1)

my$Dz_I_P=$U[8]-$U[$j];			# +value: DONOR side: distance from PREVIOUS exon end (Ez)
my$D1_I_P=$U[8]-($U[$i+1]+1);	# -value: ACCEPTOR side: distance to NEXT exon start (E1)

# E1 (first), E2 (second), Ey (second last), Ez (last) positions on exon:
# my $E_1_P = $U[$i] + 1;# my $E_2_P = $U[$i] + 2;# my $E_y_P = $U[$j] - 1;# my $E_z_P = $U[$j];

# Iy (second last), Iz (last), I1 (first) I2 (second) position on intron:
# my $I_1_P = $U[$j] + 1;# my $I_2_P = $U[$j] + 2;# my $I_y_P = $U[$i + 1] - 1;# my $I_z_P = $U[$i + 1];

		if(($U[8]>$U[$i])&&($U[8]<($U[$j]+1))&&$U[24]=~/\+/){print"$Az_E_P\t$D1_E_P\tE\t$Exon_P\t$_"}		
		if(($U[8]>$U[$j])&&($U[8]<($U[$i+1]+1))&&$U[24]=~/\+/){print"$Dz_I_P\t$D1_I_P\tI\t$Exon_P\t$_"}

# ----------------- MINUS STRAND --------------------------
# $i          $j        $i+1       $j+1    
# EzEy...E2E1<-IzIy...I2I1<-EzEy...E2E1

# MINUS DIRECTION <------
# ------ Exon ------|______ INTRON______|------ Exon ------
#                Acceptor             Donor

my$Az_E_M=$U[8]-($U[$j]+1);	# -value: ACCEPTOR side: distance from PREVIOUS intron end (Iz)
my$D1_E_M=$U[8]-$U[$i];		# +value: DONOR side: distance to NEXT intron start (I1)

my$Dz_I_M=$U[8]-($U[$i]+1);	# -value: DONOR side: distance from PREVIOUS exon end (Ez)
my$A1_I_M=$U[8]-$U[$j-1];	# +value: distance to NEXT exon start (E1)

# E1 (first), E2 (second), Ey (second last), Ez (last) positions on exon:
# my $E_1_M = $U[$j];# my $E_2_M = $U[$j] - 1;# my $E_y_M = $U[$i] + 2;# my $E_z_M = $U[$i] + 1;

# Iy (second last), Iz (last), I1 (first) I2 (second) position on intron:
# my $I_1_M = $U[$i];# my $I_2_M = $U[$i] - 1;# my $I_y_M = $U[$j - 1] + 2;# my $I_z_M = $U[$j - 1] + 1;

		if((($U[8]>$U[$i])&&($U[8]<($U[$j]+1)))&&$U[24]=~/\-/){print"$Az_E_M\t$D1_E_M\tE\t$Exon_M\t$_"}
		if((($U[8]>$U[$j-1])&&($U[8]<($U[$i]+1)))&&$U[24]=~/\-/){print"$Dz_I_M\t$A1_I_M\tI\t$Exon_M\t$_"}}}close(IN);close(O)}

# -------------- STEP2: Determine the variations that cross splice sites -------------
			open(IN,"$CAT_R/A\_$FILE");
				while(<IN>){my@G=split/\t|\,/;

					my$R_L=length($G[14]);		# [] - calculated in upstream workflow
					# my $A_L = length($G[15]); # [] - calculated in upstream workflow
					
					# unlike exon/intron starts and ends (see below figure)...
					# $From and $To are never 0; they are either greater or smaller than ZERO:
					my$from=abs($G[0]); 		# distance FROM PREVIOUS intron/exon
					my$to=abs($G[1]);   		# distance TO NEXT intron/exon
					# my $nmer = $G[24];		# 250-mer element position
					
					# my $Left_Ref = substr($G[14],0,1);
					# my $Left_Alt = substr($G[15],0,1);
					my$s=$G[28];
					
					# length of intron/exon
					my$inex_L=abs($to)+abs($from)-1;	# VALIDATE THIS!

# ------------ ACROSS SPLICE JUNCTION VARIANTS: they may or may not affect splicing and/or protein sequence -----------------------------
						# Define CROSS-JUNCTION VARIATIONS:
						# After translations, the plus strand intronic variants affecting exons must be corrected by +1 on affected exons

						# Intron and Exon CROSS-JUNCTION variants
						my$X_J=((($R_L>$to)&&$s=~/\+/&&$R_L>=2)||(($R_L>$from)&&$s=~/\-/&&$R_L>=2));				
						# Intron variants that mutate (II2|Iyz)
						my$I_var_eq=(($R_L==$to)&&$s=~/\+/&&$G[2]eq'I')||(($R_L==$from)&&$s=~/\-/&&$G[2]eq'I');
						my$I_var_1=((($R_L+1)==$to)&&$s=~/\+/&&$G[2]eq'I')||((($R_L+1)==$from)&&$s=~/\-/&&$G[2]eq'I');

						# ALL CODING and NON-CODING CROSS-JUNCTION VARIANTS ARE INCLUDED
						if($X_J){open(JUNCTION,">>$CAT_R/JIE\_$FILE");select(JUNCTION);		# INTRON (I) and EXON (E) splice junction (J) variants (JIE)
							print;close(JUNCTION)}											# if Exon variant then goes into TRANSLATION
																							# intron variants, may go into translation -> see the next step						
						# Below two files alter the conserved bases on intron (I12yz)
						elsif($I_var_eq){open(Ieq,">>$CAT_R/I\_1z\_$FILE");select(Ieq);		# INTON variants reaching or at INTRON start I1 or end Iz 
							print "$inex_L\t\t\t\t\t$_";close(Ieq)}							# they don't cross splice junction	
																							# they DON'T go into TRANSLATION; will be REPORTED
																							# NEGATIVE TRANSLATION CONTROL: does not produce protein level results
	
						elsif($I_var_1){open(Iun,">>$CAT_R/I\_2y\_$FILE");select(Iun);		# INTRON variants reaching or at INTRON start second base (I2) or end second last base (Iy)
							print "$inex_L\t\t\t\t\t$_";close(Iun)}							# they don't cross splice junction or affect I1 or Iz
																							# they DON'T go into TRANSLATION; will be REPORTED
																							# NEGATIVE TRANSLATION CONTROL: does not produce protein level results
							
						# Variants that do not cross splice junction and do not alter INTRON structure
						# Coding Exon variants go into TRANSLATION
						elsif(!$X_J&&!$I_var_eq&&!$I_var_1){open(O,">>$CAT_R/b\_$FILE");select(O);
						print "$inex_L\t\t\t\t\t$_";close(O)}}close(IN);

# STEP3: determine the variations that actually alter GT..AG structure or alternative I12..Iyz from the J_$FILE
# LEFT-match for all variants is either 1 or 0.
# If (Ref_L - LEFT_match) > ($From|$To) then (I12|Iyz) is mutated

# make additional output files:
my@out_f=('JIE_pix','JIE_pex','JIE_pni','JIE_pne','JIE_mix','JIE_mex','JIE_mni','JIE_mne');
foreach my$out_f(@out_f){open(out_f,">>$CAT_R/$out_f\_$FILE");close(out_f)}

			open(IN,"$CAT_R/JIE\_$FILE");
				while(<IN>){my@U=split/\t|\,/; # first comma denotes the start of UCSC GB data
# ---------------- PLUS STRAND ----------------------------------------
# EXON:   +value: ACCEPTOR side: distance from PREVIOUS intron end (Iz)
# EXON:   -value: DONOR side: distance to NEXT intron start (I1)
# INTRON: +value: DONOR side: distance from PREVIOUS exon end (Ez)
# INTRON  -value: ACCEPTOR side: distance to NEXT exon start (E1)
# 
#  PLUS
#  ------->
#  - - - intron - - - ->|- - - - - - EXON - - - - - ->|- - - - intron - - - ->|- - - - - EXON - - - -
#  . . . . . . . . . . .|1 2 3 . vE. . . . . . .-3-2-1|1 2 3 . . vI . . -1-2-3| . . . . . . . . . . .
#  . . . . . . . . . . .|$From . . . . . . . . . . $To|$From . . . . . . . $To| . . . . . . . . . . .
# 
# -------------- INTRON --------------------
# ACCEPTOR SIDE: =MID(V1,(98 - C1),2) where:
# C1 is $To and V1 is 250mer, and (98 of Excel = 97 of Perl) -> formula determines Iyz
# DONOR SIDE (NOT NEEDED for JUNCTION-CROSSING variants): =MID(V1,101-B1,4) where: -> formula determines I12
# B1 is $From and V1 is 250mer and 101 of (Excel = 100 of Perl)
# -------------- EXON ---------------------
# DONOR SIDE: =MID(V1,(98 - C1),2) where:
# C1 is $To and V1 is 250mer, and (98 of Excel = 97 of Perl) -> formula determines Eyz-> for I12, add +2
# ACCEPOR SIDE (Not needed for JUNCTION-CROSSING variants): =MID(V1,101-B1,4) where: -> formula determines E12 -> for Iyz, add +2
# -------------------------------------------------------------------------------------------------------------------------------

					my$Nmer=$U[24];
					my$From=$U[0];my$To=$U[1];
					my$Ref=$U[14];my$Ref_L=length($U[14]);
					my$Alt=$U[15];my$Alt_L=length($U[15]);
					# my $DIF = $Ref_L - $Alt_L;
					my$Cod=$U[2];
					my$S=$U[28];

					my$INEX_L=abs($To)+abs($From)-1;

					# $J in any variable stands for JUNCTION
					my$J_O_P=substr($Nmer,(97-$To-2),7);		# If INTRON then IwxyzE123; Ref, Length is 7 -> alter if needed
																# If EXON then  EwxyzI123; Ref, Length is 7 -> alter if needed

					my$J_O_M=substr($Nmer,(97-$From-2),7);		# If INTRON then I4321Ezyx; Ref, Length is 7 -> alter if needed
																# If EXON then E4321Izyx; Ref, lenght is 7 -> alter if needed

					# Joint: $O_front$REF_in_Nmer$O_back MUST MATCH the original 250-mer: qc pass 01/31/2016
					# '$O_' stands for the original Nmer
					my$O_front=substr($Nmer,0,99);
					my$REF_in_Nmer=substr($Nmer,99,$Ref_L);		# MUST EXACTLY MATCH w/ $Ref: qc pass 01/31/2016
					my$O_back=substr($Nmer,(99 + $Ref_L),-1);

					# Don't use perl join-function to reconstruct Nmer -> produces wrong results
					# Simple concatenation as shown below does it correctly
					my$RECONSTR_Ref_Nmer="$O_front$REF_in_Nmer$O_back";		# use as a control if needed
					my$RECONSTR_Alt_Nmer="$O_front$Alt$O_back";				# use as a control if needed

					my$J_V_P=substr($RECONSTR_Alt_Nmer,(97 - $To - 2),7);	# If INTRON then IwxyzE123; Alt, Length is 7 -> alter if needed
																			# If EXON then EwxyzI123; Alt, Length is 7 -> alter if needed
							
					my$J_V_M=substr($RECONSTR_Alt_Nmer,(97 - $From - 2),7);	# If INTRON then I4321Ezyx; Alt, Length is 7 -> alter if needed
																			# If EXON then E4321Izyx; Ref, lenght is 7 -> alter if needed

					my$AG_O_P=substr($J_O_P,2,2);my$AG_V_P = substr($J_V_P,2,2);							
					my$AG_O_M=substr($J_O_M,2,2);my$AG_V_M = substr($J_V_M,2,2);

					if(($Cod eq'I')&&($S=~/\+/)&&($AG_O_P ne $AG_V_P)){open(O,">>$CAT_R/JIE_pix\_$FILE");select(O);
					print "$INEX_L\t$AG_O_P\t$AG_V_P\t$J_O_P\t$J_V_P\t$_";close(O)}
					elsif(($Cod eq'E')&&($S=~/\+/)&&($AG_O_P ne $AG_V_P)){open(O,">>$CAT_R/JIE_pex\_$FILE");select(O);
					print "$INEX_L\t$AG_O_P\t$AG_V_P\t$J_O_P\t$J_V_P\t$_";close(O)}
					elsif(($Cod eq'I')&&($S=~/\+/)&&($AG_O_P eq $AG_V_P)){open(O,">>$CAT_R/JIE_pni\_$FILE");select(O);
					print "$INEX_L\t$AG_O_P\t$AG_V_P\t$J_O_P\t$J_V_P\t$_";close(O)}
					elsif(($Cod eq'E')&&($S=~/\+/)&&($AG_O_P eq $AG_V_P)){open(O,">>$CAT_R/JIE_pne\_$FILE");select(O);
					print "$INEX_L\t$AG_O_P\t$AG_V_P\t$J_O_P\t$J_V_P\t$_";close(O)}
					elsif(($Cod eq'I')&&($S=~/\-/)&&($AG_O_M ne $AG_V_M)){open(O,">>$CAT_R/JIE_mix\_$FILE");select(O);
					print "$INEX_L\t$AG_O_M\t$AG_V_M\t$J_O_M\t$J_V_M\t$_";close(O)}
					elsif(($Cod eq'E')&&($S=~/\-/)&&($AG_O_M ne $AG_V_M)){open(O,">>$CAT_R/JIE_mex\_$FILE");select(O);
					print "$INEX_L\t$AG_O_M\t$AG_V_M\t$J_O_M\t$J_V_M\t$_";close(O)}
					elsif(($Cod eq'I')&&($S=~/\-/)&&($AG_O_M eq $AG_V_M)){open(O,">>$CAT_R/JIE_mni\_$FILE");select(O);
					print "$INEX_L\t$AG_O_M\t$AG_V_M\t$J_O_M\t$J_V_M\t$_";close(O)}
					elsif(($Cod eq'E')&&($S=~/\-/)&&($AG_O_M eq $AG_V_M)){open(O,">>$CAT_R/JIE_mne\_$FILE");select(O);
					print "$INEX_L\t$AG_O_M\t$AG_V_M\t$J_O_M\t$J_V_M\t$_";close(O)}}
					close(IN);unlink("$M_R/$FILE");unlink("$CAT_R/A\_$FILE");
close (D_m);close(D_cat);


# ----------------------- Merge --------------------------------------
opendir(D_cat,"$CAT_R");opendir(D_a,"$A_R/$FILE");opendir(D_m,"$M_R");

open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
print "\nSPLICE\_$FILE - all entries; must be equal to $FILE\n";
print "b\_$FILE - conserved intron residues NOT affected; variant does not cross splice junction\n";
print "I\_1z\_$FILE and I\_2y\_$FILE - variant ends in first\/last and second first\/last base on intron respectively\n";
print "JIE\_$FILE plus I\_1z\_$FILE plus I\_2y\_$FILE; must add up to SPLICE\_$FILE\n";
print "JIE\_$FILE - variants across splice junction; may mutate conserved intron residues or modify or reconstitute splice-junctions\n";
print "JIE\_prefix\_$FILE - different flavors of JIE\_$FILE; must add up to JIE\_$FILE\n\n";close(O);

open(IN,"$CAT_R/JIE\_$FILE");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
while(<IN>){}print "JIE\_$FILE,$.\n";close(IN);close(Z);unlink("$CAT_R/JIE\_$FILE");

my@file=grep{/^\w/}readdir D_cat;
foreach my$file(@file){open(IN,"$CAT_R/$file");open(ALL,">>$CAT_R/SPLICE\_$FILE");select(ALL);
while(<IN>){print "$file\t$_"}close(IN);close(ALL);copy("$CAT_R/SPLICE\_$FILE","$A_R/$FILE/SPLICE\_$FILE");copy("$CAT_R/SPLICE\_$FILE","$M_R/$FILE")}
close(D_cat);close(D_a);close(D_m);
###################################
# Summarize and delete:
opendir(D_cat,"$CAT_R");opendir(D_a,"$A_R/$FILE");
my@file=grep{/^\w/}readdir D_cat;
foreach my$file(@file){open(IN,"$CAT_R/$file");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
while(<IN>){}print "$file,$.\n";close(IN);close(Z);unlink("$CAT_R/$file")}
close(D_cat);close(D_a);

# my $finis = <STDIN>;