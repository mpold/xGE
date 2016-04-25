#!/usr/bin/perl -w
use Data::Dumper;
use POSIX;
use File::Copy;
use Parallel::ForkManager;

my$pm=new Parallel::ForkManager(4);

my$MASTER_R='';
my$T_R='';
my$TR_R='';
my$A_R='';
my$M_R='';

foreach(my$NYX=1;$NYX<=24;$NYX++){
opendir(D_temp,"$T_R");
open(O,">>$T_R/$NYX\.pl");select(O);
print"#!/usr/bin/perl -w\nuse Data::Dumper;\nuse POSIX;\n";
print"my\$T_R='$T_R';\nmy\$TR_R='$TR_R';\n";
print"my\$Tx_R=\"\$TR_R\/v\_$NYX\";\nmy\$C_R=\"\$TR_R\/v\_$NYX\/Comparison\";\n";
print"my\$A_R=\"$A_R\";\n";
print"my\$M_R=\"$M_R\";\n";

open(I,"$MASTER_R/z_VARIANT_TRANSLATE.pl");while(<I>){print if!/^(\#|\t{1,}\#)/}close(I);close(O);
close(D_temp)}

# execute reference translation in parallel (see above)
# if files are very large (thousands of entries) then the below parallel execution does not help - not enough memory in 8-core PC
opendir(D_temp,"$T_R");
	my@file=grep{/\.pl$/}readdir D_temp;
		foreach my$file(@file){
$pm->start and next;
		system("$T_R/$file");
		unlink("$T_R/$file");
$pm->finish;
		}
$pm->wait_all_children;
close(D_temp);



# print"\nPost-translation and -alignent chromsome totals\n";
# print"They should match chromosome totals before translation (see above)\n";