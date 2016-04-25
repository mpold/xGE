#!/usr/bin/perl -w
use Data::Dumper;
use POSIX;
use File::Copy;
use Parallel::ForkManager;

my$pm=new Parallel::ForkManager(8);

# --- PARALLELIZE ----------
foreach(my$i=1;$i<=24;$i++){
$pm->start and next;

my$C_R="///v_$i/Comparison";

# ---------------------------------
opendir(D_c,"$C_R");
	my@file=grep{/^\d/}readdir D_c;
		foreach my$file(@file){
			open(I,"$C_R/$file");
				while(<I>){
					if($.==1){@L_1=split/ /;
						$E=@L_1;				# number of elements in array (each line as an array has a different number of elements)
						$J=$L_1[3];				# element that determines splice junction involvement
						$from=$L_1[9];
						$to=$L_1[10];										
						$EI=$L_1[11];			# E or I: special cases needs to be taken into account (see previous scripts)				
						$q=$L_1[12];			# variant exon/intron	
						$Ref=$L_1[23]; 			# reference sequence
						$Alt=$L_1[24]; 			# alternate sequence
						$Ref_L=length($Ref);
						$Alt_L=length($Alt);
						$S=$L_1[37];			# strand
						$Tx_S=$L_1[38]; 		# TxStart
						$Tx_E=$L_1[39]; 		# TxEnd
						$C_St=$L_1[40]; 		# cdsStart
						$C_E=$L_1[41]; 			# cdsEnd
						$Ec=$L_1[42]; 			# Exon total
						$L_Dif=$Alt_L-$Ref_L;	# positive when insertion; negative when deletion; 0 when non-length variant
						$Cx_Tx_E=$Tx_E+$L_Dif;
						$Cx_CDS_E=$C_E+$L_Dif;			
				}}close(I);

	#------------ CORRECT EXOME COORDINATES -----------------------
	# PLUS strand: $PEXNE refers to the 'pex' and 'pne' prefixes from JUNCTION.pl (see the script for details)
	# MINUS strand: $MEXNE refers to the 'mex' and 'mne' prefixes from JUNCTION.pl (see the script for details)
	# PLUS strand: $PIXNI refers to the 'pix' and 'pni' prefixes from JUNCTION.pl (see the script for details)
	# MINUS strand: $MIXNI refers to the 'mix' and 'mni' prefixes from JUNCTION.pl (see the script for details) 
	$PEXNE=(($EI eq'E')&&($S=~/\+/)&&($Ec>1)&&($J=~/JIE\_/));
	$MEXNE=(($EI eq'E')&&($S=~/\-/)&&($Ec>1)&&($J=~/JIE\_/));
	$MIXNI=(($EI eq'I')&&($S=~/\-/)&&($Ec>1)&&($J=~/JIE\_/));
	$PIXNI=(($EI eq'I')&&($S=~/\+/)&&($Ec>1)&&($J=~/JIE\_/));
	
	# -- no length variants other than splice junction MNV ----
	$No_L=$L_Dif==0&&!($J=~/JIE\_/);
	# -- plus and minus strand length variation ---
	$PLUS=($S=~/\+/&&!($J=~/JIE\_/)&&($L_Dif ne 0));
	$MINUS=($S=~/\-/&&!($J=~/JIE\_/)&&($L_Dif ne 0));
	# -- singel exon gene: the same for both minus and plus strand ----
	$SINGLE=($Ec==1&&($L_Dif ne 0));

if($No_L){}	# do nothing if a no-length, no-splice junction MNV variation
# ---- common for all splice junction variants ----
if($PEXNE||$MEXNE||$MIXNI||$PIXNI||$PLUS ||$MINUS){
	# Do LINE #1: split ends and starts into separate lines
	open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
		while(<I>){if($.==1){
			print"@L_1[43..(43+$Ec-1)]\,";					# exon starts
			print"\n";
			print"@L_1[(43+$Ec+1)..(43+(2*$Ec))]\,";		# exon ends
			print"\n";
	}}close(I);close(O);
}
# ------ PLUS strand non-splice junction length variants --------
if($PLUS){
	if($Ec>1){
			open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
				while(<I>){
				if($.==4){@ps_4=split/ |\,/;
					print"@ps_4[0..($q-1)]\,"; 						# exons starts not corrected
						# ----- corrected exon starts ------
						foreach(my$a=($q);$a<=($Ec-1);$a++){
						my$CxP_s=$ps_4[$a]+$L_Dif;print"$CxP_s\,"} 	# exon starts corrected by $L_Dif (rest of the exons)
						print"\,\n";
					}
				elsif($.==5){@ps_5=split/ |\,/;
					print"@ps_5[0..($q-2)]\,"; 						# exon ends not corrected (one less than starts)
						# -- corrected exon ends -----------
						foreach(my$b=($q-1);$b<=($Ec-1);$b++){
						my$CxP_e=$ps_5[$b]+$L_Dif;print"$CxP_e\,"} 	# exon ends corrected by $L_Dif (rest of the exons)
						print"\,\n";						
			}}close(I);close(O);
		}}
# ------ MINUS strand non-splice junction length variants -------
if($MINUS){		
	if($Ec>1){
			open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
				while(<I>){
				# -- the whole array is [0..($Ec-1)]: both starts and ends
				if($.==4){@ms_4=split/ |\,/;
					$UX=($Ec-1)-($q-1); 		# upper limit for uncorrected; note that uncorrected starts are -1
												# - keep $UX as a GLOB
					print"@ms_4[0..($UX)] ";	# uncorrected ends
					# -------- corrected ends -------------						
					foreach(my$c=($UX+1);$c<=($Ec-1);$c++){
					my$CXM_E=$ms_4[$c]+$L_Dif;print"$CXM_E\,"} 	# exon ends corrected by $L_Dif
					print"\,\n";
					}
				elsif($.==5){@ms_5=split/ |\,/;
					print"@ms_5[0..($UX-1)] ";					# uncorrected starts
					# -------- corrected start -------------------------------------						
					foreach(my$d=($UX);$d<=($Ec-1);$d++){
					my$CXM_S=$ms_5[$d]+$L_Dif;print"$CXM_S\,"} 	# exon starts corrected by $L_Dif
					print"\,\n";
			}}close(I);close(O);
		}}
# ----------- add annotation line and format ----------------	
if($PLUS||$MINUS){
		# --------- reconstitute annotations ----------------
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
				if($.==1){
					print"@L_1[0..37] ";
					print"$Tx_S ";			# Tx_start
					print"$Cx_Tx_E ";		# Tx_end
					print"$C_St ";			# Cds start
					print"$Cx_CDS_E ";		# Cds end
					print"$Ec ";			# Exon count
					print"\n";
					}
				if($.==6){s/^( |\,)//;s/( \,)$/\,/;s/(\, )$/\,/;print}
				if($.==7){s/^( |\,)//;s/( \,)$/\,/;s/(\, )$/\,/;print}
				# --------- reconstitute top line -------------------
				if(eof){print"@L_1[(43+(2*$Ec)+2)..((43+(3*$Ec)+6))]"}
		}close(I);close(O);
		# ------------- final format --------------------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);while(<I>){
		if($.>7&&$.<12){s/\,/ /g;chomp;print}}print"\n";close(I);close(O);
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);
		while(<I>){if($.==2){print}}print"\n";close(I);close(O);
		unlink("$C_R/$file");rename("$C_R/P\_$file","$C_R/$file");
		}		
# --------- all single exon genes, strad independent -----------
if($SINGLE){
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);
			while(<I>){
				if($.==1){my$S_E_P=$L_1[45]+$L_Dif;		# single exon gene exon end
					print"@L_1[0..38] ";				# left side of Line #1; includes Tx_start
					print"$Cx_Tx_E ";					# corrected Tx_end
					print"$L_1[40] ";					# CDS start
					print"$Cx_CDS_E ";					# corrected CDS end
					print"$Ec ";						# exon count
					print"$L_1[43]  ";					# first exon start
					print"$S_E_P ";						# corrected first exon end
					print"@L_1[(44+(2*$Ec))..($E-2)] ";	# rest of Line #1
					print"\n";
					}
			elsif($.==2){print;print"\n"}
		}close(I);close(O);
		unlink("$C_R/$file");rename("$C_R/P\_$file","$C_R/$file");
	}

# ------ SPICE JUNCTION CONDIONS ----------------------------
if($PEXNE){
		# -- the whole array is [0..($Ec-1)]: starts and ends
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
			if($.==4){@p_4=split/ |\,/;
				print"@p_4[0..($q-1)]\,"; 						# exons starts not corrected
					# ----- corrected exon starts ------
					foreach(my$m=($q);$m<=($Ec-1);$m++){my$Cx_s=$p_4[$m]+$L_Dif;print"$Cx_s\,"} 	# exon starts corrected by $L_Dif (rest of the exons)
					print"\,\n";
				}
			elsif($.==5){@p_5=split/ |\,/;
				print"@p_5[0..($q-2)]\,"; 						# exon ends not corrected (one less than starts)
				my$f_p_e=$p_5[($q-1)]+($to+1);					# first corrected exon end, by $to+1
				print"$f_p_e\,";
					# -- the rest of modified exon ends ---
					foreach(my$n=($q);$n<=($Ec-1);$n++){my$Cx_e=$p_5[$n]+$L_Dif;print"$Cx_e\,"} 	# exon ends corrected by $L_Dif (rest of the exons)
					print"\,\n";						
		}}close(I);close(O);	
	}

if($MEXNE){
		# -------- corrections -----------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
			# -- the whole array is [0..($Ec-1)]: both starts and ends
			if($.==4){@m_4=split/ |\,/;
				$Ux=($Ec-1)-($q-1); 	# upper limit for uncorrected; note that uncorrected starts are -1
										# - keep $Ux as a GLOB
				print"@m_4[0..$Ux] ";	# uncorrected ends
				# -------- corrected ends -------------						
				foreach(my$z=($Ux+1);$z<=($Ec-1);$z++){
				my$CX_E=$m_4[$z]+$L_Dif;print"$CX_E\,"} 	# exon ends corrected by $L_Dif
				print"\,\n";
				}
			elsif($.==5){@m_5=split/ |\,/;
				print"@m_5[0..($Ux-1)] ";					# uncorrected starts
				# -------- corrected start -------------------------------------						
				my$f_P_s=$m_5[($Ux)]+($from+1);				# first corrected exon start, by $from+1
				print"$f_P_s\,";

				foreach(my$zz=($Ux+1);$zz<=($Ec-1);$zz++){my$CX_S=$m_5[$zz]+$L_Dif;print"$CX_S\,"} 	# exon starts corrected by $L_Dif
				print"\,\n";
		}}close(I);close(O);
	}
if($PEXNE||$MEXNE){
		# --------- reconstitute annotations ----------------
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
				if($.==1){
					print"@L_1[0..37] ";
					print"$Tx_S ";			# Tx_start
					print"$Cx_Tx_E ";		# Tx_end
					print"$C_St ";			# Cds start
					print"$Cx_CDS_E ";		# Cds end
					print"$Ec ";			# Exon count
					print"\n";
					}
				if($.==6){s/^( |\,)//;print}
				if($.==7){s/^( |\,)//;print}
				# --------- reconstitute top line -------------------
				if(eof){print"@L_1[(43+(2*$Ec)+2)..((43+(3*$Ec)+6))]"}
		}close(I);close(O);
		# ------------- final format --------------------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);while(<I>){
		if($.>7&&$.<12){s/\,/ /g;chomp;print}}print"\n";close(I);close(O);
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);
		while(<I>){if($.==2){print}}print"\n";close(I);close(O);
		unlink("$C_R/$file");rename("$C_R/P\_$file","$C_R/$file");
		}
if($MIXNI){
		# -- the whole array is [0..($Ec-1)]: both starts and ends
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
			my$un=($Ec-($q));
			# --- correct exon ends ---
			if($.==4){@M_4=split/ |\,/;
				print"@M_4[0..($un-1)] ";					# uncorrected exon ends				
				#--- FIRST end correction ---
				my$f_m_e=$M_4[$un]+($from+1);
				print"$f_m_e\,";
				# ---- rest of corrected ends ---------
				foreach(my$x=($un+1);$x<=($Ec-1);$x++){my$Cx_E=$M_4[$x]+$L_Dif;print"$Cx_E\,"} 	
				print"\,\n";
				}
			# --- correct exon starts ----
			elsif($.==5){@M_5=split/ |\,/;
				print"@M_5[0..($un-1)] ";					# uncorrected exon starts
				# -- FIRST start correction------
				my$f_m_s=$M_5[$un]+($L_Dif);
				print"$f_m_s\,";
				# ---- rest of corrected starts -------	
				foreach(my$y=($un+1);$y<=($Ec-1);$y++){my$Cx_S=$M_5[$y]+$L_Dif;print"$Cx_S\,"}
				print"\,\n";
			}}close(I);close(O);
		# -- reconstitute the top line annotations, w/ corrected coordinates ----
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);while(<I>){
		# --- use only IF, no ELSIF ---------	
			if($.==1){
				print"@L_1[0..37] ";
				print"$Tx_S ";			# Tx_start
				print"$Cx_Tx_E ";		# Tx_end
				print"$C_St ";			# Cds start
				print"$Cx_CDS_E ";		# Cds end
				print"$Ec ";			# Exon count
				print"\n";
				}
			if($.==6){s/^(\,| )//;print}
			if($.==7){s/^(\,| )//;print}
			# --------- fully reconstitute top line -------------
			if(eof){print"@L_1[(43+(2*$Ec)+2)..((43+(3*$Ec)+6))]"}
		}close(I);close(O);	
		# ------------- final format --------------------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);while(<I>){
		if($.>7&&$.<12){s/\,/ /g;chomp;print}}print"\n";close(I);close(O);
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);
		while(<I>){if($.==2){print}}print"\n";close(I);close(O);
		unlink("$C_R/$file");rename("$C_R/P\_$file","$C_R/$file");		
	}
if($PIXNI){
		# -- the whole array is [0..($Ec-1)]: both starts and ends
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
			while(<I>){
				# ----- exon starts -------
				if($.==4){@P_4=split/ |\,/;
					print"@P_4[0..($q-1)] "; 			# exons starts not corrected
					my$F_s_p=$P_4[$q]+($to+1);			# FIRST corrected exon start
					print"$F_s_p\,";
												
					# -- rest of corrected exon starts ---
					foreach(my$e=($q+1);$e<=($Ec-1);$e++){my$cx_s=$P_4[$e]+$L_Dif;print"$cx_s\,"}
					print"\,\n";
					}
				# ----- exon ends ------------
				elsif($.==5){@P_5=split/ |\,/;
					print"@P_5[0..($q-1)] ";			# uncorrected exon ends
					my$F_e_p=$P_5[$q]+$L_Dif;			# FIRST corrected exon end
					print"$F_e_p\,";	
					# -- rest of corrected exon ends -----
					foreach(my$s=($q+1);$s<=($Ec-1);$s++){my$cx_e=$P_5[$s]+$L_Dif;print"$cx_e\,"}
					print"\,\n";
		}}close(I);close(O);

		# -- reconstitute annotations, w/ corrected coordinates -----
		open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);while(<I>){
		# --------- use only IF, no ELSIF ---------	
			if($.==1){
				my$Cx_EI='e';			# used if PIXNI, adjustment from Intron to Exon
				my$Cx_q=$q+1;			# used if PIXNI, adjustment for exon number
				print"@L_1[0..10] ";
				print"$Cx_EI ";
				print"$Cx_q ";
				print"@L_1[13..37] ";
				print"$Tx_S ";
				print"$Cx_Tx_E ";
				print"$C_St ";
				print"$Cx_CDS_E ";
				print"$Ec ";
				print"\n";
				}
			if($.==6){s/^(\,| )//;print}
			if($.==7){s/^(\,| )//;print}
			# --------- fully reconstitute top line -------------
			if(eof){print"@L_1[(43+(2*$Ec)+2)..((43+(3*$Ec)+6))]"}
		}close(I);close(O);	
		# ------------- final format --------------------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);while(<I>){
		if($.>7&&$.<12){s/\,/ /g;chomp;print}}print"\n";close(I);close(O);
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);
		while(<I>){if($.==2){print}}print"\n";close(I);close(O);
		unlink("$C_R/$file");rename("$C_R/P\_$file","$C_R/$file");
		}		
	}
close(D_c);
$pm->finish;
}
$pm->wait_all_children;
