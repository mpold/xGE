# OUTPUT: coding and non-coding variants separated
#
# ARCHIVED FILES:
# 
# 1. C_$FILE - contains all variants that go into translation
# 2. nr_C_$FILE - the non-redundant version of the above file. They should be the same though
# 3. NON_$FILE - contains all variants that do not go into translation
# 4. nr_NON_$FILE - the non-redundant version of NON_$FILE. Both files should be the same
# 
# no elements added to the array from the previous step, which as PolyA (see 7_x_POLYA.pl)


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
	my @file=grep{/^(\w|\d)/}readdir D_m;
		foreach my $file(@file){$FILE=$file;
			open(IN,"$M_R/$file");
				while(<IN>){my@CODING = split /\t|\,/;
# Elements come from the array from the previous step (PolyA)
# [2] - SPL = JIE (any JIE containing acronym)
# [10] - EI_id = E
# [14] - TAG_1 = atg_u5pC, stop_u3mC, ncRNA_P, ncRNA_M
# [15] - MAP = CDS_P,CDS_M
					my$SPL=$CODING[2]=~/JIE\_/&&$CODING[15]=~/CDS\_/;
					my$TAG_1=$CODING[14]=~/atg\_u5pC/||$CODING[14]=~/stop\_u3mC/;
					my$MAP=$CODING[15]=~/CDS\_/&&$CODING[10]eq'E';

					# print variants that go into translation into the same file
					# if needed they can be separated
					if($SPL){open(C,">>$CAT_R/C\_$FILE");select(C);print;close(C)}
					elsif($TAG_1){open(C,">>$CAT_R/C\_$FILE");select(C);print;close(C)}
					elsif($MAP){open(C,">>$CAT_R/C\_$FILE");select(C);print;close(C)}
					# print non-coding variants in a separate file
					elsif(!$SPL && !$TAG_1 && !$MAP){open(NON,">>$CAT_R/NON\_$FILE");select(NON);print;close(NON)}}
			close(IN)}
close(D_m);close(D_c);
# ------------------ CONTROL NR --------------------------------------------------------------------------------------
# Output files from the previous step should both be non-redundant but just in case, a control nr-step is taken here..
# ..to make sure they indeed are so

opendir(D_c,"$CAT_R");
	my@file=grep{/\w/}readdir D_c;
	foreach my $file(@file){open(NR,"$CAT_R/$file");open(nr,">>$CAT_R/nr\_$file");select(nr);
	my %seen;while(<NR>){print if$seen{$_}++==0}close(NR);close(nr)}
close(D_c);

# -------- SUMMARIZE ------------------------------- 
opendir(D_c,"$CAT_R");
open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
print"\nC\_$FILE - all entries that go into translation\n";
print"NON\_$FILE - non\-coding entries; not to be translated\n";
print"nr\_C\_$FILE should equal C\_$FILE; and nr\_NON\_$FILE should equal NON\_$FILE\n";
close(O);

my@file=grep{/\w/}readdir D_c;
	foreach my$file(@file){open(IN,"$CAT_R/$file");open(O,">>$A_R/$FILE/Zummary\_$FILE\.csv");select(O);
	while(<IN>){}print"$file,$.\n";close(IN);close(O);copy("$CAT_R/$file","$A_R/$FILE/$file")}
close(D_c);

opendir(D_c,"$CAT_R");
	my@file_x=grep{!/^C\_/}readdir D_c;
	foreach my$file_x(@file_x){unlink("$CAT_R/$file_x");rename("$CAT_R/C\_$FILE","$CAT_R/$FILE")}
close(D_c);