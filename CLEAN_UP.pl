#!/usr/bin/perl -w
use 5.020;
use strict;
use warnings;

use Data::Dumper;
use POSIX;
use Parallel::ForkManager;
my $pm = new Parallel::ForkManager(4);

my $Temp_R = '';
my $M_R = '';

my $R_R = '';
my $T_R = '';
my $CAT_R = '';

# EMPTY FOLDERS:
# Enter the directory variables as an array
my @Directory = ("$M_R","$R_R","$CAT_R","$Temp_R");
	foreach my $Directory (@Directory){
		opendir (DIR,"$Directory");
			my @file = grep {/./} readdir DIR;
				foreach my $file (@file){
			unlink ("$Directory/$file");
			}
	close(DIR);
}

# remove DNA motif files and folder:
	opendir (DIR,"$Temp_R/$Chr");			
		my @file = grep {/./} readdir DIR;
			foreach my $file (@file){
				unlink ("$Temp_R/$Chr/$file");		
			}
	close(DIR);
rmdir("$Temp_R/$Chr");

# remove parallelized indel folders sub-folders in Temporary root ($Temp_R):
# upper $i corresponds to the number of parallels in SNV pipeline
foreach (my $i = 1; $i <= 64; $i++){
	opendir (DIR,"$Temp_R/$i");			
		my @file = grep {/./} readdir DIR;
			foreach my $file (@file){
				unlink ("$Temp_R/$i/$file");		
			}
	close(DIR);
rmdir("$Temp_R/$i");
}

# remove folders and sub-folders in $T_R:
# upper $i corresponds to the number of parallels in pipeline
foreach (my $i = 1; $i <= 64; $i++){
	opendir (DIR,"$T_R/v\_$i/$Comparison");			
		my @file = grep {/./} readdir DIR;
			foreach my $file (@file){
				unlink ("$T_R/v\_$i/$Comparison/$file");		
			}
	close(DIR);
rmdir("$T_R/v\_$i/$Comparison");

# remove SNV folders:
	opendir (DIR,"$T_R/v\_$i");			
		my @file = grep {/./} readdir DIR;
			foreach my $file (@file){
				unlink ("$T_R/v\_$i/$file");		
			}
	close(DIR);
rmdir("$T_R/v\_$i");

# remove translation pre-curson files
	opendir (DIR,"$T_R");			
		my @file = grep {/./} readdir DIR;
			foreach my $file (@file){
				unlink ("$T_R/$file");		
			}
	close(DIR);
}
