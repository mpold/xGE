# FUTURE: integrate polyA in UTR3
# Good point to start: http://dnafsminer.bic.nus.edu.sg
# 
# OUTPUT FILES: 
# Prefixes: self-explanatory (see Zummary_$FILE)
#
# ARCHIVED FILES: goes into next step
# 1. PA\_$FILE: contains both polyA consensus matching and non-matching entries:
##### APPLIES TO file 1. and 2. PolyA-matching entries contain values in the added tabs: ref- and alt-polyA consensus (see the body of this script for additional detail)
# 2. PolyA_$FILE: contains entries that contain (potentially) variant polyA consensus
# 3. No_polyA_$FILE: contains entries that do not contain variant polyA consensus
#
# OUTPUT ELEMENTS: 33 + UCSC exome descriptions
# Elements from this script:
# [0] - pA_REF: contains $i_polyA_consensus of reference sequence (polyA_Ref); value = '' if no match with polyA consensus found
# [1] - pA_ALT: mutated polyA_Ref; value = '' if no match with polyA consensus found
# Elements from previous script:
# [2] - [32] + UCSC elements
#
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

opendir(D_m,"$M_R");opendir(D_c,"$CAT_R");
	my@file=grep{/^(\w|\d)/}readdir D_m;
		foreach my $file(@file){$FILE=$file;
			open(IN,"$M_R/$file");
				while(<IN>){
# no warnings;
					my@U=split/\t|\,/; # first comma denotes the start of UCSC GB data

					my$UTR3=$U[13];
					my$Nmer=$U[30];
					my$Exon=$U[8];
					my$PolyA_window=substr($Nmer,94,11);	# a 12-nucleotide regions, each position can be an Alt start position
					my$SIG_P='(AATAAA|ATTAAA)';				# PLUS strand poly A consensus signals
					my$SIG_M='(TTTATT|TTTAAT)';				# MINUS strand poly A consensus signals
					my$Ref=$U[20];my$Ref_L=length($Ref);
					my$Alt=$U[21];my$Alt_L=length($Alt);
					my$Str=$U[34];
	
					foreach(my$i=0;$i<=5;$i++){
						# define the length and pattern-matching area for reference polyA signal candidates
						# 11 bases, each of which can be the variant start coordinate is searched for polyA consensus (see above)
						my$Consensus_Ref=substr($PolyA_window,$i,6);
						my$Utr_3_P=($Consensus_Ref=~/$SIG_P/&&$UTR3=~/(U3\_P)/&&$Exon eq'E')if$Str=~/\+/;
						my$Utr_3_M=($Consensus_Ref=~/$SIG_M/&&$UTR3=~/(U3\_M)/&&$Exon eq'E')if$Str=~/\-/;

						if($Utr_3_P||$Utr_3_M){open(P,">>$CAT_R/PolyA\_$file");select(P);
						my$Ref_Cx=99-(94+($i-1));
						my$Ref_Front=substr($Nmer,(94+$i),($Ref_Cx-1));
						my$Ref_Rear=substr($Nmer,(94+$i+$Ref_Cx),$i);

						# Alt-polyA -- $Var_Poly_A -- is stitched together in three parts:
						# 1. $Ref_Front
						# 2. $Alt
						# 3. $Ref_Rear
						my$Var_Poly_A="$Ref_Front$Alt$Ref_Rear";
						my$POLY_A_PATTERN=($Var_Poly_A=~/$SIG_P/&&$Str=~/\+/)||($Var_Poly_A=~/$SIG_M/&&$Str=~/\-/);

						# $i in a printout indicates the start of polyA pattern found in the 12-base search window
						# On $Nmer - $i = 0 means position 94, and 5 means position 99
						# 94 means that the polyA end and variant start coincide
						# 99 means that the polyA start and variant start coincide
						print"$i\_$Consensus_Ref\t$Var_Poly_A\t$_"if!$POLY_A_PATTERN}close(P)}}close(IN);

		# Define the polyA-matching entries -> via chr.pos.ID.ref.alt
		open(P,"$CAT_R/PolyA\_$file");open(String,">>$CAT_R/query\_$file");select(String);
		while(<P>){my@poly=split/\t|\,/;print"\t$poly[19]\t$poly[20]\t$poly[21]\t$poly[22]\t$poly[23]\t\n"}
		print"nothing\n";close(String);close(P);
		# nr the above entries
		open(String,"$CAT_R/query\_$file");open(nr,">>$CAT_R/nr\_$file");select(nr);
		my%seen;while(<String>){print if$seen{$_}++==0;}close(String);close(nr);
		# make a search scalar
		open(nr,"$CAT_R/nr\_$file");open(OUT,">>$CAT_R/my\_$file");select(OUT);
		while(<nr>){s/$/\|/;chomp;print}close(nr);close(OUT);

		# finalize the search pattern in OR (polyA matching hits)
		open(OUT,"$CAT_R/my\_$file");open(O,">>$CAT_R/o\_$file");select(O);
		while(<OUT>){s/\|nothing\|//;s/\|$//;
			my$pattern=$_;
			# define as GLOBAL variable
			$PATTERN="\($pattern\)"if$.==1;
			print"$PATTERN";}close(OUT);close(O);
			# print entries that don't match polyA-pattern (in a file where poly_A variations exist):
			if($PATTERN=~/\|/){open(IN,"$M_R/$FILE");open(N,">>$CAT_R/No\_polyA\_$FILE");select(N);
			while(<IN>){print"\t\t$_" if!/$PATTERN/}close(IN);close(N);}			
			# print out everything if no polyA-pattern (in a file where poly_A variations DO NOT exist):
			if($PATTERN=~/nothing/){open(IN,"$M_R/$FILE");open(N,">>$CAT_R/No\_polyA\_$FILE");select(N);
			while(<IN>){print"\t\t$_"}close(IN);close(N)}
		unlink("$M_R/$FILE")}

# ------------- Summarize, copy and delete ----------
open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
print"\nPolyA\_$FILE: lists variants possible affecting polyadenylation signal in 3\'UTR\nYES_polyA_$FILE and NO_polyA\_$FILE must add up to PA\_$FILE\n";
close(Z);
	# count entries with variant polyA-site:
	open(IN,"$CAT_R/PolyA\_$FILE");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
	while(<IN>){}print"YES\_polyA\_$FILE,$.\n";close(IN);close(Z);
	copy("$CAT_R/PolyA\_$FILE","$A_R/$FILE/PolyA\_$FILE");
	# count entries other than variant polyA-site:
	open(IN,"$CAT_R/No_polyA_$FILE");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
	while(<IN>){}print"NO\_polyA\_$FILE,$.\n";close(IN);close(Z);
	copy("$CAT_R/No_polyA_$FILE","$A_R/$FILE/No_polyA_$FILE");
	# Merge the two above input files: MERGER will be the output into the next step:
	open(IN,"$CAT_R/PolyA\_$FILE");open(M,">>$A_R/$FILE/PA\_$FILE");select(M);
	while(<IN>){print};close(IN);
	open(IN,"$CAT_R/NO\_polyA\_$FILE");open(M,">>$A_R/$FILE/PA\_$FILE");select(M);
	while(<IN>){print};close(IN);close(M);
	# count entries in the merger:
	open(IN,"$A_R/$FILE/PA\_$FILE");open(Z,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(Z);
	while(<IN>){}print"PA\_$FILE,$.\n";close(IN);close(Z);
	
close(D_m);close(D_c);
# --------------------------------------------------
opendir(D_c,"$CAT_R");my@file=grep{/^\w/}readdir D_c;
foreach my $file(@file){unlink("$CAT_R/$file")}
copy("$A_R/$FILE/PA\_$FILE","$M_R/$FILE");
close(D_c);