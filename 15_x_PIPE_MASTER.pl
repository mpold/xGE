#!/usr/bin/perl -w
use strict;
use warnings;

use Data::Dumper;
use File::Copy;
use POSIX;


my$MASTER_R='';

system("$MASTER_R/CLEAN_UP.pl");
system("$MASTER_R/1_x_0_M1.pl");
system("$MASTER_R/2_x_adj_NMER.pl");
system("$MASTER_R/3_x_SORT_QC.pl");
system("$MASTER_R/4_x_INPUT.pl");
system("$MASTER_R/5_x_CAT_1.pl");
system("$MASTER_R/6_x_ATG_EXT.pl");
system("$MASTER_R/7_x_JUNCTION.pl");
system("$MASTER_R/8_x_POLYA.pl");
system("$MASTER_R/9_x_CODING.pl");
system("$MASTER_R/10_x_T_MASTER.pl");
system("$MASTER_R/11_x_PRE_TRANS.pl");
system("$MASTER_R/12_x_REF_TRANS.pl");
system("$MASTER_R/13_x_Cx_COORDINATE.pl");
system("$MASTER_R/14_x_ALT_TRANS.pl");

