# ------- write a perl batch-script headers: this script fetches specified exons from full-lenght transcripts ------------------
open(B,">>$Tx_R/BATCH\.pl");select(B);print"\#\!\/usr\/bin\/perl \-w\nuse Data::Dumper;\n";print"my\$Tx_R=\"$Tx_R\";\n";close(B);


# ------- REFERENCE PROTEIN TRANSLATION -----------
# Calculate exon size (downstream offset)
# Calculate the transcript and exon start positions
opendir(D_tx,"$Tx_R");
	my@file=grep{/^\d/}readdir D_tx;
		foreach my$file(@file){
			# Determine the size and position of of exons
			open(I,"$Tx_R/$file");open(O,">>$Tx_R/$file");select(O);print"\n";
				while(<I>){if($.==1){my@E=split/ /;
# [20] - chromosome
# [21] - variant position
# [33] - Nmer
# [37] - Strand
# [38] - TxStart
# [39] - TxEnd
# [40] - cdsStart
# [41] - cdsEnd
# [42] - Exon count

# ---------  Transcript and coding region coordinates ------------------------
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
								my$L=$E[$j]-$E[$e];
								# Tx start coordinate -> $E[10]
								my$E_Pos=$E[$e]-$E[38];
								# Print out PLUS strand exons, starting with the first coding
								print"$E_Pos,$L,$Up_P,\+,\n" if$E[37]=~/\+/&& 
								# Specifies the coding exons as non-negative reading frames found in UCSC annotations
								($E[$e+(2*abs($E[42]))+6]ne-1);
								# Required to determine the first coding exon
								my$M_Cx=$M_TxL-$E_Pos;my$M_CDS_CX=$M_Cx-$Up_M;
								# Print out MINUS strand exons, starting with the first coding
								print"$E_Pos,$L,$M_CDS_CX,\-,\n" if$E[37]=~/\-/&& 
								# Specifies the coding exons as non-negative reading frames found in UCSC annotations
								($E[$e+(2*abs($E[42]))+6]ne-1)}}
							}
						close(I);close(O);

# --- CORRECT FOR THE CDS START SITE AND PRINT EXON START POSTION and SIZE from GENE SEQUENCE ---------
			open(I,"$Tx_R/$file");open(O,">>$Tx_R/$file");select(O);
				while(<I>){my@First=split/ |\,/;
					my$Atg=$First[1]-abs($First[2]-$First[0]) if$.==5&&/\+\,$/;
					print"$First[2],$Atg,cds,\+\n" if$.==5&&/\+\,$/;
					print"$First[0],$First[1],cds,\+\n" if$.>5&&/\+\,$/;
					# Define the last line in the file as 'eof'
					my$Last=$_ if eof;
					print"$First[0],$First[1],cds,\-\n"if$.>4&&/\-\,$/&&!$Last;
					print"$First[0],$First[2],cds,\-\n"if$Last&&/\-\,$/;
				}close(I);close(O);

# --- write temporary perl file to print out exons starting with the ATG of the first coding exon -------
# the 3'-s is printed out entirely
# don't write a separate perl file for each file in a directory - becomes way too slow
			open(I,"$Tx_R/$file");open(O,">>$Tx_R/BATCH\.pl");select(O);
				print"open(I,\"\$Tx_R/$file\");open(O,\">>\$Tx_R/$file\");select(O);";
				print"while(<I>){if(\$.==2){";

				while(<I>){my@Start_End=split/\,/ if$.>5;
				if(/,cds,(\+|\-)$/&&$.>5){print"my\$M\_$.=substr(\$_,$Start_End[0],$Start_End[1]);print\"\$M\_$.\";"}}
			print"}}close(I);close(O);";close(I);close(O)}
close(D_tx);

# Execute each chromosome-specific transcript-assembly script:
system("$Tx_R/BATCH\.pl");
unlink("$Tx_R/BATCH\.pl");

# ---------- REFERENCE PROTEIN TRANSLATION --------
# Calculate exon size (downstream offset)
# Calculate the transcript and exon start positions
opendir(D_tx,"$Tx_R");
	my@file=grep{/^\d/}readdir D_tx;
		foreach my$file(@file){
			# Determine the size and position of of exons
			open(I,"$Tx_R/$file");open(O,">>$Tx_R/$file");select(O);print"\n";close(O);
# --------- FORMAT THE EXONS FOR TRANSLATION ----------------------------------------------------	
			open(O,">>$Tx_R/1\_$file");select(O);
			while(<I>){chomp if$.>5&&!/\,/&&!/^$/;print$_ if$.>5&&!/\,/&&!/^$/}close(I);close(O);
			
# -------- Input file content -------------------------------------------------------------------------
# 1\_$file (TX-file handle) - only the assembled transcript (a single line without a newline at the end)
# $file (I-filehandle) - (1) variant description, (2) full Tx (unspliced), (3) all coding exon size data, (4) spliced ATG-TxEnd
# 2\_$file (O-filehandle) - formatted and protected Alt-sequence (MINUS only) spliced ATG-TxEnd

			open(TX,"$Tx_R/1\_$file");open(O,">>$Tx_R/2\_$file");open(I,"$Tx_R/$file");select(O);
				while(<I>){
	
					if($.==1){my@Strand=split/ /;
						# GLOB - strand: PLUS(+) or MINUS(-)
						$S=$Strand[37]}
					
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
			unlink("$Tx_R/1\_$file");

			# FIND FIRST STOP
			open(I,"$Tx_R/2\_$file");open(O,">>$Tx_R/3\_$file");select(O);
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
			unlink("$Tx_R/2\_$file");				

			# print the ATG-STOP into a separate file:
			open(I,"$Tx_R/3\_$file");open(O,">>$Tx_R/4\_$file");select(O);
			while(<I>){print if$.==1}close(I);close(O);
			unlink("$Tx_R/3\_$file");				

			open(I,"$Tx_R/4\_$file");open(O,">>$Tx_R/P\_$file");select(O);
				while(<I>){
					# Determine the codon count:
					# last number denotes the stop codon
					my$N=0;
					$N+=scalar(split(/\,/));
					# print out codon numbers:
					foreach(my$i=1;$i<=$N;$i++){print"$i "}print"\n";
					my@Cod=split/\,/;
					
# --------------- TRANSLATE REFERENCE TRANSCRIPTS --------------------
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
unlink("$Tx_R/4\_$file");
		# -------- Format ---------------------------------------------------
		open(I,"$Tx_R/$file");open(O,">>$Tx_R/P\_$file");select(O);print"\n";
		while(<I>){print if!/^$/;}close(I);close(O);
		unlink("$Tx_R/$file");rename("$Tx_R/P\_$file","$Tx_R/$file");
}
close(D_tx);


