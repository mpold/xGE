# OUTPUT of this script  that goes into NEXT STEP: N_$file -> $FILE
# Tab is used as a delimiter then: the same as in the previous step
# [0] CHROM:	chromosome - e.g. chr1
# [1] CX_POS:	corrected position - e.g. 26609354, corrected if match found at 0-position only (see [5])
# [2] ID:		period (.) or dbSNP_ID (comes from .vcf)
# [3] REF:		Reference sequence, e.g. CCT
# [4] ALT:		Alternative sequence, e.g. C
# [5] MATCH:	0 or 1 or both match -> 1, 2, or 21	
# [6] TYPE:		variation type - S (SNV), M (MNV), or L (length)
# [7] ORI_POS:	original position from the .vcf file - 26609354
# [8] DOM:		if starts with N, then a 0-to-1-based adjustment carried out (to unify downstream, NMER-based decisions)
#				# adding an N simply "pushes" the 0-based position "downstream" by one base

# ARCHIVED FILES: 
# 1. N_$file (N_$FILE)
# 2. 1_$FILE:				all original zero or both 0-and-1-based entries
# 3. 0_$FILE:				all original 0-based entries
# 4. Zummary_$FILE.csv:		line counts of N_$FILE, 1_$FILE, and 0_$FILE

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

# --------- QC, and variant inclusion and exclusion --------------------------------------
opendir(D_i,"$I_R");
	my @file=grep{/^\w/}readdir D_i;
		foreach my$file(@file){$FILE=$file;
			# make the rest of the files in $I_R for the current foreach:
			open(ONE,">>$I_R/1\_$file");close(ONE);open(ZERO,">>$I_R/0\_$file");close(ZERO);
				open(IN,"$I_R/$file");
					while(<IN>){my@VCF=split/\t/;my$P=$VCF[5];
						if(($P eq 2)||($P eq 21)){open(ONE,">>$I_R/1\_$file");select(ONE);print;close(ONE)}
						elsif($P eq 1){open(ZERO,">>$I_R/0\_$file");select(ZERO);
						s/$VCF[8]/N$VCF[8]/;	# add N to each 0-based NMER
						print;close(ZERO)}}close(IN);unlink("$I_R/$file");
}
close(D_i);
# --------- SUMMARIZE AND ARCHIVE ----------------------
opendir(D_i,"$I_R");opendir (D_a,"$A_R/$FILE");
open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
print"\n1- and 0-based totals; they must add up to N\_$FILE\n";close(OUT);
	@file=grep{/^\d/}readdir D_i;
		foreach my $file(@file){open(IN,"$I_R/$file");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
		while(<IN>){}print "$file,$.\n";close(IN);close(OUT);		
		copy("$I_R/$file","$A_R/$FILE/$file");
		# --------------- MERGE ----------------------------------
		open(IN,"$I_R/$file");open(OUT,">>$I_R/$FILE");select(OUT);
		while(<IN>){print}close(IN);close(OUT);
		copy("$I_R/$FILE","$A_R/$FILE/N\_$FILE");unlink("$I_R/$file");
		}
		open(IN,"$A_R/$FILE/N\_$FILE");open(OUT,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(OUT);
		while(<IN>){}print "N\_$FILE,$.\n";close(IN);close(OUT);
close(D_i);close(D_a);
