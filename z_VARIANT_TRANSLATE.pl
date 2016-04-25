# OUTPUT:
# LINE 1: reference transcript description including UCSC transcript coordinates
# LINE 2: variant transcript description including UCSC transcript coordinates in which length variant coordinates are adjusted
# LINE 3: ref transcript amino acid/codon numbers
# LINE 4: ref transcript amino acids
# LINE 5: ref transcript codons
# LINE 6: variant transcript amino acid/codon numbers
# LINE 7: variant transcript amino acids
# LINE 8: variant transcript codons


# write perl-batch-file header:
open(B,">>$C_R/BATCH\.pl");select(B);print"\#\!\/usr\/bin\/perl \-w\nuse Data::Dumper;\n";print"my\$C_R=\"$C_R\";\n";close(B);

#######################################
##### VARIANT PROTEIN TRANSLATION #####
# Calculate exon size (downstream offset)
# Calculate the transcript and exon start positions
opendir(D_c,"$C_R");
	my@file=grep{/^\d/}readdir D_c;
		foreach my$file(@file){
			# Determine the size and position of of exons
			open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);print"\n";
				while(<I>){
	
				if($.==1){my@E=split/ /;
# [20] - chromosome
# [21] - variant position
# [33] - Nmer
# [37] - Strand
# [38] - TxStart
# [39] - TxEnd
# [40] - cdsStart
# [41] - cdsEnd
# [42] - Exon count
### Transcript and coding region coordinates ##############
# Distance from transcript start to CDS start (PLUS strand)
my$Up_P=$E[40]-$E[38];
# Distance from TxStart to CDS start (MINUS strand)
my$Up_M=$E[39]-$E[41];
# the entire length of MINUS stand transcript
my$M_TxL=$E[39]-$E[38];
						# number of exons (coding and non-coding) = $E[42]
						# $e is a constant - element [] of the first exon start
						foreach(my$e=43;$e<=($E[42]+42);$e++){
							my$j=$e+($E[42]+1);
								my$Length=$E[$j]-$E[$e];
								# Tx start coordinate -> $E[10]
								my$E_Pos=$E[$e]-$E[38];
								# Print out PLUS strand exons, starting with the first coding
								print"$E_Pos,$Length,$Up_P,\+,\n" if$E[37]=~/\+/&& 
								# Specifies the coding exons as non-negative reading frames found in UCSC annotations
								($E[$e+(2*abs($E[42]))+6]ne-1);
								# Required to determine the first coding exon
								my$M_Cx=$M_TxL-$E_Pos;my$M_CDS_CX=$M_Cx-$Up_M;
								# Print out MINUS strand exons, starting with the first coding
								print"$E_Pos,$Length,$M_CDS_CX,\-,\n" if$E[37]=~/\-/&& 
								# Specifies the coding exons as non-negative reading frames found in UCSC annotations
								($E[$e+(2*abs($E[42]))+6]ne-1)}}}
						close(I);close(O);

# --------- CORRECT FOR THE CDS START SITE AND PRINT EXON START POSTION and SIZE from GENE SEQUENCE ------------------
			open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);
				while(<I>){my@First=split/ |\,/;
					my$Atg=$First[1]-abs($First[2]-$First[0]) if$.==5&&/\+\,$/;
					print"$First[2],$Atg,cds,\+\n" if$.==5&&/\+\,$/;
					print"$First[0],$First[1],cds,\+\n" if$.>5&&/\+\,$/;
					# Define the last line in the file as 'eof'
					my$Last=$_ if eof;
					print"$First[0],$First[1],cds,\-\n"if$.>4&&/\-\,$/&&!$Last;
					print"$First[0],$First[2],cds,\-\n"if$Last&&/\-\,$/;
				}close(I);close(O);

# ----- write temporary perl file to print out exons starting with the ATG of the first coding exon -------------------
# the 3'-s is printed out entirely
# don't write a separate perl file for each file in a directory - becomes way too slow
			open(I,"$C_R/$file");open(O,">>$C_R/BATCH\.pl");select(O);
				print"open(I,\"\$C_R/$file\");open(O,\">>\$C_R/$file\");select(O);";
				print"while(<I>){if(\$.==2){";

				while(<I>){my@Start_End=split/\,/ if$.>5;
				if(/,cds,(\+|\-)$/&&$.>5){print"my\$M\_$.=substr(\$_,$Start_End[0],$Start_End[1]);print\"\$M\_$.\";"}}
			print"}}close(I);close(O);";close(I);close(O)}
close(D_c);
# ------- Execute each chromosome-specific transcript-assembly script ---------
system("$C_R/BATCH\.pl");
unlink("$C_R/BATCH\.pl");

# --------- VARIANT PROTEIN TRANSLATION -----------
# Calculate exon size (downstream offset)
# Calculate the transcript and exon start positions
opendir(D_c,"$C_R");
	my@file=grep{/^\d/}readdir D_c;
		foreach my$file(@file){
			# Determine the size and position of of exons
			open(I,"$C_R/$file");open(O,">>$C_R/$file");select(O);print"\n";close(O);

# --------- FORMAT THE EXONS FOR TRANSLATION ------	
			open(O,">>$C_R/1\_$file");select(O);
			while(<I>){chomp if$.>5&&!/\,/&&!/^$/;print$_ if$.>5&&!/\,/&&!/^$/}close(I);close(O);

# Input file content:
# 1\_$file (TX-file handle) - only the assembled transcript (a single line without a newline at the end)
# $file (I-filehandle) - (1) variant description, (2) full Tx (unspliced), (3) all coding exon size data, (4) spliced ATG-TxEnd
# 2\_$file (O-filehandle) - formatted and protected Alt-sequence (MINUS only) spliced ATG-TxEnd

			open(TX,"$C_R/1\_$file");open(O,">>$C_R/2\_$file");open(I,"$C_R/$file");select(O);
				while(<I>){
	
					if($.==1){
						my@Strand=split/ /;
						# GLOB - strand: PLUS(+) or MINUS(-)
						$S=$Strand[37];
						}
					
					while(<TX>){
						# add characters as needed
						# DO THE PLUS STRAND FORMATTING: order of substitutions not relevant:
						if($S=~/\+/){s/T/T /g;s/t/t /g;s/A/A /g;s/a/a /g;s/C/C /g;s/c/c /g;s/G/G /g;s/g/g /g;
							s/X/X /g;s/x/x /g;s/N/N /g;s/n/n /g;s/U/U /g;s/u/u /g; # X,N,U may not be necessary
						print}

						if($S=~/\-/){
						# DO THE MINUS STRAND FORMATTING
						# ORDER of SUBSTITIONS MUST BE FOLLOWED AS BELOW:
						# a - o
						# c - p
						# g - q
						# t - r
						# STEP1: Protect the lower case bases (Alt-sequence) first:
						s/a/ o/g;s/c/ p/g;s/g/ q/g;s/t/ r/g;					
						# STEP2: Complement the rest of the sequence:
						s/T/ a/g;s/A/ t/g;s/C/ g/g;s/G/ c/g;s/X/ x/g;s/N/ n/g;s/U/ a/g;
						# STEP3: remake the large caps
						s/t/T/g;s/a/A/g;s/c/C/g;s/g/G/g;s/x/X/g;s/n/N/g;s/u/U/g;
						# STEP3: unprotect, and complement the Alt-sequence
						s/o/t/g;s/p/g/g;s/q/c/g;s/r/a/g;
						# STEP4: reverse minus strand transcripts
				print scalar reverse}}}print "\n";close(TX);close(O);close(I);
			unlink("$C_R/1\_$file");

			# FIND FIRST STOP
			open(I,"$C_R/2\_$file");open(O,">>$C_R/3\_$file");select(O);
				while(<I>){my$C=0;
					# every third white space is replaced by comma
					s/ /(++$C%3==0)?",":$&/ge;s/ //g;
					# Newline very first stop codon
					# ATG-to-STOP placed in a separate line
					# since stop-codons can be a mixed upper-lower case motifs (stop-codon variants)..
					#.. the subs must be done one all possible upper-lowe case combinations:
s/TAG/TAG\n/;s/tag/tag\n/;s/Tag/Tag\n/;s/tAg/tAg\n/;s/taG/taG\n/;s/TAg/TAg\n/;s/tAG/tAG\n/;s/TaG/TaG\n/;
s/TGA/TGA\n/;s/tga/tga\n/;s/Tga/Tga\n/;s/tGa/tGa\n/;s/tgA/tgA\n/;s/TGa/TGa\n/;s/tGA/tGA\n/;s/TgA/TgA\n/;
s/TAA/TAA\n/;s/taa/taa\n/;s/Taa/Taa\n/;s/tAa/tAa\n/;s/taA/taA\n/;s/TAa/TAa\n/;s/tAA/tAA\n/;s/TaA/TaA\n/;	
			print}close(I);close(O);
			unlink("$C_R/2\_$file");				

			# print the ATG-STOP into a separate file:
			open(I,"$C_R/3\_$file");open(O,">>$C_R/4\_$file");select(O);
			while(<I>){print if$.==1}close(I);close(O);
			unlink("$C_R/3\_$file");				

			open(I,"$C_R/4\_$file");open(O,">>$C_R/P\_$file");select(O);
				while(<I>){
					# Determine the codon count:
					# last number denotes the stop codon
					my$N=0;
					$N+=scalar(split(/\,/));
						# print out codon numbers:
						foreach(my$i=1;$i<=$N;$i++){print"$i "}print"\n";
						my@Cod=split/\,/;
					
# ---------- TRANSLATE ALT TRANSCRIPTS ---------------------------------
					foreach my$Cod(@Cod){

print"F "if$Cod=~/(TTT|TTC)/i;
print"L "if$Cod=~/(TTA|TTG|CTT|CTC|CTA|CTG)/i;
print"I "if$Cod=~/(ATT|ATC|ATA)/i;
print"M "if$Cod=~/(ATG)/i;
print"V "if$Cod=~/(GTT|GTC|GTA|GTG)/i;
print"S "if$Cod=~/(TCT|TCC|TCA|TCG|AGT|AGC)/i;
print"P "if$Cod=~/(CCT|CCC|CCA|CCG)/i;
print"T "if$Cod=~/(ACT|ACC|ACA|ACG)/i;
print"A "if$Cod=~/(GCT|GCC|GCA|GCG)/i;
print"Y "if$Cod=~/(TAT|TAC)/i;
print"X "if$Cod=~/(TAA|TAG|TGA)/i;
print"H "if$Cod=~/(CAT|CAC)/i;
print"Q "if$Cod=~/(CAA|CAG)/i;
print"N "if$Cod=~/(AAT|AAC)/i;
print"K "if$Cod=~/(AAA|AAG)/i;
print"D "if$Cod=~/(GAT|GAC)/i;
print"E "if$Cod=~/(GAA|GAG)/i;
print"C "if$Cod=~/(TGT|TGC)/i;
print"W "if$Cod=~/(TGG)/i;
print"R "if$Cod=~/(CGT|CGC|CGA|CGG|AGA|AGG)/i;
print"G "if$Cod=~/(GGT|GGC|GGA|GGG)/i;
			}print"\n@Cod"}close(I);close(O);
unlink("$C_R/4\_$file");
		# ------- Finish all required formats ----------------------------
		open(I,"$C_R/$file");open(O,">>$C_R/P\_$file");select(O);print"\n";
		while(<I>){print if!/^$/;}close(I);close(O);
		unlink("$C_R/$file");
		
		rename("$C_R/P\_$file","$C_R/$file");
}
close(D_c);
# ------------ END OF VARIANT NM TRANSLATION ----------

# ------------ MERGE, TRANSLATE & ALIGN ---------------
# MERGE the translated wild-type and variant transcript
opendir(D_c,"$C_R");opendir(D_tx,"$Tx_R");my@file=grep{/^\d/}readdir D_tx;

	# to count all Ref vs. Alt alignment files at the end of this loop (see the very bottom)
	my$count=0;

foreach my$file(@file){# print reference UCSC description
open(W,"$Tx_R/$file");open(M,">>$Tx_R/M\_$file");select(M);
while(<W>){print if$.==5}close(W); # close(M);
# print variant description - coordinates of the length variants modified
open(V,"$C_R/$file");open(M,">>$Tx_R/M\_$file");select(M);
while(<V>){print if$.==5}close(V); # close(M);
# print reference translation; turn the plus strand large cap into small cap; add a space at the end of codons
open(W,"$Tx_R/$file");open(M,">>$Tx_R/M\_$file");select(M);
while(<W>){s/$/ / if$.==3;print if$.==1||$.==2||$.==3}close(W); # close(M);	
# print variant translation; turn the plus strand large cap into small cap; add a space at the end of codons
open(V,"$C_R/$file");open(M,">>$Tx_R/M\_$file");select(M);
while(<V>){s/$/ / if$.==3;print if$.==1||$.==2||$.==3}print"\n";close(V);close(M);
# clean up
unlink("$Tx_R/$file");unlink("$C_R/$file");rmdir("$C_R");rename("$Tx_R/M\_$file","$Tx_R/$file");

# ---------- START OF ALIGNMENT -------------------------------------------------------------------------------
# Determine the number of elements in the aligment lines in order to adjust the downstream steps appropriately:
open(I,"$Tx_R/$file");open(O,">>$Tx_R/N\_$file");select(O);while(<I>){no warnings;
# see the top of this script for the line descriptions
if($.==3||$.==4||$.==5||$.==6||$.==7||$.==8){my@STACK=split/ /;my$Nu=$#STACK+1;print"$Nu,"}}close(I);close(O);
# temporarily switch to non-descriptive alignment lines only
open(I,"$Tx_R/$file");open(O,">>$Tx_R/A\_$file");select(O);while(<I>){no warnings;
# see the top of this script for the line descriptions
if($.==3||$.==4||$.==5||$.==6||$.==7||$.==8){print}}close(I);close(O);

# --------- VERY IMPORTANT STEP WITHOUT WHICH THE REST WON'T WORK PROPERLY ----------------
# OUTPUT has a format that enables all lines be of the same length for the downstream steps	
open(N,"$Tx_R/N\_$file");open(A,"$Tx_R/A\_$file");open(O,">>$Tx_R/E\_$file");select(O);
while(<N>){no warnings;
my@SEQ=split/\,/;
# Determine the highest and lowest number in @SEQ
my($low,$high)=(sort{$a<=>$b}@SEQ)[0,-1];
while(<A>){my@REF=split/ /;print"@REF[0..($high)]\n"}}
close(N);close(O);close(A);unlink("$Tx_R/N\_$file");unlink("$Tx_R/A\_$file");

# Format the file from the previous step into a file which all lines contain the same number of delimiters
# Delimiter = white space
open(E,"$Tx_R/E\_$file");open(K,">>$Tx_R/K\_$file");select(K);
while(<E>){chomp if/(\w|\d)/;print}print"\n";close(E);close(K);unlink("$Tx_R/E\_$file");

# CONTROL to verify that the number of elements in each scalar is the same
open(K,"$Tx_R/K\_$file");open(C,">>$Tx_R/C\_$file");select(C);while(<K>){no warnings;
my@KL=split/ /;my$M=@KL;					
# see the top of the script for line descriptions; subtract 2 because the descriptive lines are not included
print"$M,"if($.==1||$.==2||$.==3||$.==4||$.==5||$.==6)}close(K);close(C);

# Format K_$file file from the previous step into a scalar 
# Delimiter = white space		
open(K,"$Tx_R/K\_$file");open(S,">>$Tx_R/S\_$file");select(S);
while(<K>){chomp;print}close(K);close(S);unlink("$Tx_R/K\_$file");

# builds a pair-wise comparison portion in the original merger file
open(S,"$Tx_R/S\_$file");open(C,"$Tx_R/C\_$file");open(O,">>$Tx_R/$file");select(O);	
	while(<C>){my@Num=split/\,/;
	
	# control condition make sure the upstream equalization of alignment elements works, otherwise no alignment
	# GLOBAL
	$HIGH=$Num[1]if$Num[0]==$Num[1]&&$Num[1]==$Num[2]&&$Num[2]==$Num[3]&&$Num[3]==$Num[4]&&$Num[4]==$Num[5];
	while(<S>){no warnings;my@MISMATCH=split/ /;
	print"alignment_OK\n\n"if$HIGH;
	print"upstream_failures\n\n"if!$HIGH
}}close(S);close(C);close(O);unlink("$Tx_R/S\_$file");unlink ("$Tx_R/C\_$file");

# count files in each chromosome specific Tx folder:
$count++if$file=~/^\d/;
}

	# obtain the input file name from the mapping folder:
	opendir(D_m,"$M_R");my@Input=grep{/^(\w|\d)/}readdir D_m;
		foreach my$Input(@Input){$FILE=$Input}	# get the original input-file name in global
	close(D_m);

	# KEEP the commented line for now - it prints out the file names chromosome/Tx_R specifically
	# open(O,">>$A_R/$FILE/z_Post_translation_Summary.csv");select(O);print"$Tx_R,@file,\n";close(O);
	open(O,">>$A_R/$FILE/Zummary\_$FILE.csv");select(O);print"$Tx_R,$count\n";close(O);

close(D_c);close(D_tx);

