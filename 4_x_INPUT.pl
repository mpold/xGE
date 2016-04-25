# x_INPUT.pl joins variants with appropriate transcripts
# ARCHIVE: m_file, contains all variant-transcript joins non-redundantly
# 
# OUTPUT: variant-NM/NR-joins, goes into next step: m_$FILE
# TOTAL of ELEMENTS: 15 + UCSC entry elements
# [0] - [14] from x_SORT.pl; all elements now tab-separated followed by..
# .. NR- or NM-record from UCSC Genome Browser download
#
# ARCHIVED FILE = V_Tx_$FILE

#!/usr/bin/perl -w
# use strict;
use warnings;

use Data::Dumper;
use POSIX;
use Parallel::ForkManager;
use File::Copy;


my$pm=new Parallel::ForkManager(8);

my$Temp_R='';
my$M_R='';
my$EXOME_R='';
my$I_R='';
my$A_R='';

# EXOME coordinate file: coordinates from UCSC Table browser
my$Tx='UCSC';

# File extensions:
my$PL = '.pl';

# Prefixes: fixed prefixes - don't change them
my$NR='NR';
my$PRE='PRE';


# ----------- STEP1: make the variant-UCSC_Tx joins ------------------
opendir(D_e,"$EXOME_R");opendir (D_i,"$I_R");opendir(D_temp,"$Temp_R");	
	my@file=grep{/^\w/}readdir D_i;
		foreach my$file(@file){$FILE=$file;
		foreach(my$i=1;$i<= 24;$i++){
$pm->start and next;
			open(IN,"$I_R/$file");
				while(<IN>){s/,X,/,23,/;s/,Y,/,24,/;
				my@Parallel=split /\,/;
				my$In_Chr=$Parallel[1];
				my$POS=$Parallel[2];

				if($In_Chr==$i){open(EX,"$EXOME_R/$i");open(O,">>$Temp_R/$i");select(O);
				while(<EX>){my@Tx=split /\t|chr/;
				print"@Parallel,$_"if($POS>$Tx[5])&&($POS<($Tx[6]+1));
				}close(EX);close(O)}}close(IN);
		$pm->finish}
$pm->wait_all_children}
close(D_temp);close(D_i);close(D_e);

# ----- STEP2: Format and merge chromosome specific variant-UCSC_Tx joins --------------
opendir(D_m,"$M_R");opendir(D_temp,"$Temp_R");opendir(D_a,"$A_R/$FILE");	
	@file=grep{/^\d/}readdir D_temp;
		foreach my$file(@file){
			open(IN,"$Temp_R/$file");open(OUT,">>$M_R/$FILE");select(OUT);
				while(<IN>){		
					s/ /\t/g if!/\,/;	# restore tabs as delimiters
					chomp if!/\,/;		# make single line entries
					print;
					}close(IN);close(OUT);unlink("$Temp_R/$file")}
	copy("$M_R/$FILE","$A_R/$FILE/V\_Tx\_$FILE");
close(D_temp);close(D_m);close(D_a);

# -------------- summarize output ----------
opendir(D_m,"$M_R");opendir(D_a,"$A_R/$FILE");	
	open(IN,"$M_R/$FILE");open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
	print"\nV\_Tx\_$FILE: all possible variant\-transcript (NM or NR) joins\n";close(O);
	open(IN,"$M_R/$FILE");open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
	while(<IN>){}print "V\_Tx\_$FILE,$.\n";close(IN);close(O);
close(D_m);close(D_a);
