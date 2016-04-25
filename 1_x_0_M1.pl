# OUTPUT of this script  that goes into NEXT STEP: M\_i\_nr\_$file -> $FILE
# Tab is used as a delimiter then:
# [0] CHROM:	chromosome - e.g. chr1
# [1] CX_POS:	corrected position - e.g. 26609354, corrected if match found at 0-position only (see [5])
# [2] ID:		period (.) or dbSNP_Id (comes from .vcf)
# [3] REF:		Reference sequence, e.g. CCT
# [4] ALT:		Alternative sequence, e.g. C
# [5] MATCH:	0 or 1 or both match, coded as 1, 2, and 21, respectively
# [6] TYPE:		variation type - S (SNV), M (MNV), or L (length)
# [7] ORI_POS:	original position from the .vcf file - 26609354
# [8] DOM:		DNA domain, 100th (-1 or both) position of which matches the first position of Ref - 250mer

# ARCHIVED FILES: 
# 1. original input file (OIF = $FILE) - just the name as it first appears in INPUT
				 ## may contain the .vcf comment lines; if so then don't use for cumulative stats

# 2. q_$FILE - good quality entries in the OIF
# 3. x1_$FILE - poor quality entries in OIF
# 4. M_i_nr_q_$FILE - contains all entries from q_OIF that match GRCh (first base of variation matches expected base in reference genome)
					  ## becomes named as OIF in INPUT when this script finishes (after OIF is archived)
# 5. x_i_nr_1_$FILE - the provided REF not matching GRCh at either 0- or 1-position; excluded from the further analysis
					  # FUTURE: design a separate work-flow for these entries
# 6. Zummary_$FILE.csv - line counts from all files described above (1.-5.)
		 
#!/usr/bin/perl -w
# use strict;
use warnings;

use Data::Dumper;
use POSIX;
use File::Copy;

my$Temp_R='';
my$CHR_R='';
my$COM_R='';
my$I_R='';
my$Out_D='';
my$A_R='';

my$Chr='chr';
my$FA='.fa';
my$PL='.pl';

mkdir("$Out_D");

# STEP1: QC and variant inclusion and exclusion
opendir(D_i,"$I_R");
	my@file=grep{/^\w/}readdir D_i;
		foreach my$file(@file){$FILE=$file;	# define input as global variable
mkdir("$A_R/$FILE");

			open(IN,"$I_R/$file");			
			# make the rest of the files in $I_R for the current foreach:
			open(BAD,">>$I_R/x1\_$file");close(BAD);
			open(GOOD,">>$I_R/q\_$file");close(GOOD);			
				while(<IN>){
					s/^chrX/chr23/;			# reverse this when trimming and qc is done!
					s/^chrY/chr24/;			# reverse this when trimming and qc is done!
					# s/^chrM/chr25/;		# FUTURE: add mitochondrial reference sequence to human genome source file
				
					my@VCF=split/chr|\t/;
					my$Ref=$VCF[4];my$Alt=$VCF[5];
					my$Ref_L=length($VCF[4]);my$Alt_L=length($VCF[5]);

#### EXCLUDE POOR QUALITY VARIANTS: 
# no difference between Ref and Alt
# chromosomes described as no 1-22, X or y
# chromosomes containing positions larger than the largest possible on Chromosome 1
# no PASS variants
# variants containing bases other than A, C, G, T, U, X, N

my$NO_VARIATION=$Ref eq $Alt; 		# they happen
my$ZERO=($Ref_L==0)||($Alt_L==0);	# Ref or Alt sequence equals blank (these will not work in TRANSLATION)

my$CHR=$VCF[1];
my$W='\D';
my$L='\W';								# if letters other than X or Y are present in original input
my$NON_NUMB=$_ if/^chr$W/;				# FUTURE: make this condition more strict
my$NON_LETT=$_ if 
($Ref =~ /$L/||$Alt=~/$L/);				# no no-letter characters in Ref or Alt
my$COMMENT=$_ if/^\#/;					# .vcf column header and comment lines

my$UN_VCF = 							# describes .vcf format requirements; if not met -> poor quality entry
($VCF[2]<250000000&&					# FUTURE: make the '<' condition chromosome specific 
$VCF[3]=~/(\.|\w)/&& 					# Current version considers only the rough length of chromosome 1 (largest chromosome)
$Ref=~/(A|C|G|T|U|N|X)+/&& 
$Alt=~/(A|C|G|T|U|N|X)+/&& 
$VCF[6]=~/\d+/&& 
$VCF[7]eq'PASS'&& 
$VCF[8]=~/AO/&& 
$CHR<25)&& 								# FUTURE: $CHR<= 25 capable to handle mitochondrial genome as well
!$NON_NUMB&& 
!$NON_LETT&& 
!$COMMENT;

	if($ZERO||$NO_VARIATION||!$UN_VCF){ # exclude entries of poor quality, or no difference between Alt and Ref
		open(BAD,">>$I_R/x1\_$file");select(BAD);print if!/^\#/;close(BAD)}
	else{open(GOOD,">>$I_R/q\_$file");select(GOOD);print;close(GOOD)}
	}close(IN);	# close the original input file
		
	# Split the good quality entries chromosome specifically in order to speed up the downstream steps
	# $GENOME[1] corresponds to chromosome numbers
	open(GOOD,"$I_R/q\_$file");while(<GOOD>){
	my@GENOME=split/chr|\t/;open(SPL,">>$Temp_R/$GENOME[1]");select(SPL);print}
	close(GOOD);close(SPL);

copy("$I_R/$file","$A_R/$FILE/$file");
copy("$I_R/x1\_$file","$A_R/$FILE/x1\_$file");
unlink("$I_R/x1\_$file");unlink("$I_R/$file");
}
close(D_i);

# ----- STEP2: make the chromosome specific .pl scripts --------
opendir(D_i,"$I_R");opendir(D_t,"$Temp_R");opendir(D_o,"$Out_D");
	@file=grep{/^(q_)/}readdir D_i;
		foreach my$file(@file){
				
		my@C=grep{/^\d/}readdir D_t; # $C - any chromosome found in an input file
			foreach my$C(@C){							
				open(IN,"$I_R/$file");open(OUT,">>$Out_D/$C$PL");select(OUT);
					print"\#\!\/usr\/bin\/perl \-w\n";
					print"use Data\:\:Dumper\;\n";
					print"open\(IN,\"$CHR_R/$Chr$C$FA\"\);";
					print"open\(OUT,\"\>\>$Out_D/$C\"\);";
					print"select\(OUT\);\n";
					print"while\(\<IN>\)\{\n";
					print"my\$Mot_L\=250;\n";	# Define the size of a larger DNA motif housing the variant					
					print"my\@P=\(";			# define chromosome position for variations
						while(<IN>){
							# split VCF input
							my@vcf=split/chr|\t/;
							print"$vcf[2]," if $vcf[1]==$C;
						}	
						print"\);\n";
						print"foreach my\$P\(\@P\)\{\n";
						print"my\$Chrom\=\$_;\n";

						# define a 250-mer wherein the variant position is found
						# if needed, expand on or trim this DNA domain, e.g. to determine the number of repeats in the variant domain
						# It may or may not, however, be so. Input can be 0, -1, +1 etc.
						print"my\$Chrom_Pos_Nmer\=\$P\-100;";
						print"my\$MOTIF_N\=substr\(\$_,\$Chrom_Pos_Nmer,\$Mot_L\);\n";
						print'print"';
						print"chr$C\t";
						print"\$P\\t\$MOTIF_N\\t\\n\";";
						print"\n\}\}";
						print"close\(IN\);close\(OUT\);\n";
				close(IN);close(OUT);
			system("$Out_D/$C$PL");
		unlink("$Out_D/$C$PL");

# ----- Join .vcf variants with their larger DNA domains (see above) ------
		open(INPUT,"$Temp_R/$C");while(<INPUT>){# 1 of 24 splits of original input file
				my@In=split/\t/;
				my$CHRM=$In[0];my $POS=$In[1];	
		open(CHROMO,"$Out_D/$C");while(<CHROMO>){
				my@MAP=split/\t/;
				my$chrm=$MAP[0];my$pos=$MAP[1];my$Nmer=$MAP[2];

				if(($chrm eq $CHRM)&&($pos eq $POS)){open(JOIN,">>$Temp_R/$file");select(JOIN);
					# Truncate .vcf entries: only essential data will enter the subsequent steps
					print"$Nmer\t@In[0..4]\t\n";	# use this line to process real .vcf files
				}}}close(JOIN);close(CHROMO);close(INPUT);
	unlink("$Temp_R/$C");unlink("$Out_D/$C");unlink("$I_R/q\_$file");
	}
# ------ STEP4: nr the output from the previous step ---------------------
		open(IN,"$Temp_R/$file");open(NR,">>$Temp_R/nr\_$file");select(NR);
		my%seen;while(<IN>){print if$seen{$_}++==0}close(IN);close(NR);

# ------ STEP5: DETERMINE 1-, 0-, or unknown-based variations ----------------------
		open (IN,"$Temp_R/nr\_$file");open(ONE,">>$Temp_R/i\_nr\_$file");select(ONE);
			while (<IN>){
				s/ /\t/g;
				my@MOTIF=split/\t/;
					my$Nmer=$MOTIF[0];
					my$Ref=$MOTIF[4];my$Alt=$MOTIF[5];
					my$R_L=length($Ref);my$A_L=length($Alt);

						# Determine if either 0- or 1-position match exists
						# output = 1 and 2 mean 0- and 1-based match
						# output = 0 means no match at either 0- or 1-position
						foreach(my $i=99;$i>=98;$i--){# $i = 9 corresponds to the 1-based scenario
							my$MATCH_1=substr($MOTIF[0],$i,$R_L);
							my$j=$i-97;
								if($MATCH_1 eq $Ref){print"$j";	 # Refs matching 0 or 1 or both are tagged 0, 1, or 01, respectively;
								}
								elsif($MATCH_1 ne $Ref){print""; # Refs not matching at 0 or 1 are not tagged
								}}
						print"\tL\t$_"if$R_L ne $A_L;
						print"\tS\t$_"if(($R_L==$A_L)&&($R_L==1));
						print"\tM\t$_"if(($R_L==$A_L)&&($R_L>1));
					}
		close(IN);close(ONE);
# ------- STEP6: Adjust variant chromosome positions if:----------
# ..variant chromosome position is only ZERO AND
		open(IN,"$Temp_R/i\_nr\_$file");
		# Make the rest of the files:
		open(A,">>$Temp_R/M\_i\_nr\_$file");close(A);open(X,">>$Temp_R/x\_i\_nr\_$file");close(X);		
			while(<IN>){
				my@ADJ=split/\t/;
					my$R=$ADJ[6];		# Ref
					my$R_adj=$ADJ[4]-1;	# Variant position MINUS ONE

						# Adjust the coordinates of length variants if (see the formula below)
						# if data come from the same .vcf file then 
						# (substr($R,0,1) ne substr($R,1,1) might not even be necessary
						# To match entire REF to 250mer use the following Excel formulas
# EXCEL CONTROL FORMULA: =IF(MID(X1,(100+J1),LEN(L1))=MID(L1,1,LEN(L1)),TRUE,FALSE)
# if only the 0-position is matching then: EXCEL CONTROL FORMULA: =IF(MID(X1,(99+J1),LEN(L1))=MID(L1,1,LEN(L1)),TRUE,FALSE)
						# note that unlike Perl, Excel uses 1-based leght system -> 100 in Excel = 99 in Perl
						if($ADJ[0]==1){open(A,">>$Temp_R/M\_i\_nr\_$file");select(A);
							print"$ADJ[3]\t$R_adj\t$ADJ[5]\t$ADJ[6]\t$ADJ[7]\t$ADJ[0]\t$ADJ[1]\t$ADJ[4]\t$ADJ[2]\t\n"}
						elsif($ADJ[0]==2||$ADJ[0]==21){open(A,">>$Temp_R/M\_i\_nr\_$file");select(A);
							print"$ADJ[3]\t$ADJ[4]\t$ADJ[5]\t$ADJ[6]\t$ADJ[7]\t$ADJ[0]\t$ADJ[1]\t$ADJ[4]\t$ADJ[2]\t\n"}
						# Exclude variants w/ Ref that does not match up with GRCh at either 0- or 1-position
						elsif($ADJ[0]==0){open(X,">>$Temp_R/x\_i\_nr\_$file");select(X); # separate non-matching Ref entries
							print"$ADJ[3]\t$ADJ[4]\t$ADJ[5]\t$ADJ[6]\t$ADJ[7]\t$ADJ[0]\t$ADJ[1]\t$ADJ[4]\t$ADJ[2]\t\n"}
						else{open(X,">>$Temp_R/x\_i\_nr\_$file");select(X); # separate non-matching Ref entries
							print"$ADJ[3]\t$ADJ[4]\t$ADJ[5]\t$ADJ[6]\t$ADJ[7]\t$ADJ[0]\t$ADJ[1]\t$ADJ[4]\t$ADJ[2]\t\n"}}
		close(IN);close(A);close(X);
copy("$I_R/$file","$A_R/$FILE/$file");
copy("$Temp_R/M\_i\_nr\_$file","$A_R/$FILE/M\_i\_nr\_$file"); # the same file is renamed in the next line and copied to the input directory
copy("$Temp_R/M\_i\_nr\_$file","$I_R/$FILE");				  # M_..file will appear in Input folder as q_file
copy("$Temp_R/x\_i\_nr\_$file","$A_R/$FILE/x\_i\_nr\_$file");
unlink("$I_R/$file"); # the 'q' prefix file
}
close(D_o);close(D_t);close(D_i);

system("$COM_R/CLEAN_UP.pl");	 # alter this as needed in the production version

# --------- STEP7: Summarize folder ---------------------------------------------
opendir(D,"$A_R/$FILE");
	@file=grep{/^\w/}readdir D;
	foreach my$file(@file){open(IN,"$A_R/$FILE/$file");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
	while(<IN>){}print "$file,$.\n" if$file ne $FILE;close(IN);close(OUT)}
	# Writes the line count of the input file at the bottom of Step 1 summary
	open(IN,"$A_R/$FILE/$FILE");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
	while(<IN>){}print "STEP\_1\: INPUT\_$FILE,$.\n\n";close(IN);close(OUT);
close(D);

