# xGE
Genome annotation pipeline for variant reporting
- The initial version is a collection of perl scripts. Each script can work separately, as well as in a pipeline.
- The script running order is in '15_x_PIPE_MASTER.pl'.
- None of the scripts in the repository contains specified directories. Hence, they need to be added. Each user can choose individual directory names.
- Currently, the scripts that write the CDS and protein level consequences following the HGSV reporting guidelines are work in progress. Will added as they mature.


SCRIPT FUNCTIONS:

1_x_0_M1:
Determines if input variants are in 0- or 1-based. All 0-based variants are adjusted to 1-based. If a variant is neither 0- nor 0-based, then it is elminated. Also, this scripts fetches a larger, 250-mer DNA domain (Nmer) within which the variation is specified in the input .vcf files. Having such a DNA domain in the pipeline for each variant will facilitate a number of steps described downstream.

2_x_adj_NMER: 
Adjusts the Nmer of all 0-based variations. A letter 'N' is added in front of each 0-based 250mer

3_x_SORT_QC: 
All non-SNV variations are QC-ed. Often, an MNV or a length variation comes described redundantly (eg. Ref: AAG to AG, where the left side A is redundant in both Ref and Alt). This type of redundancy is eliminated in by this script. At the end of this script, the genome level variation type is assigned to each variant: SNV, MNV, deletion, insertion, inversion, duplication, and complex (delins).

4_x_INPUT:
All input variants are matched to trascript coordinates obtained from UCSC Genome Browser. The following options are checked in UCSC downloads:
OPTIONS CHECKED:Genome: Human; Reference genome: GRCh37/hg19, Feb. 2009; Group: mRNA and EST; Track: Human mRNAs; Table: RefSeq Genes (refGene); Output format: all fields from the selected table

5_x_CAT_1:
Determines the downstream path for all variants: either a coding or non-coding region variant. All coding region variants will go into translation. Non-coding region variants will have their own workflows.
Additinally, this scripts evaluates if the length variations upstream of first ATG or downstream of STOP end up in the coding region. If yes, then they go into downstream translation.

6_x_ATG_EXT:
Determines if variants affecting first ATG produce inframe upstream extensions (?extM1). Also, figures out if an upstream inframe stop codon for the first ATG might exist. Combining these finding, it is possible to more accurately determine the impact of variants affecting the consensus translation start site.

7_x_JUNCTION:
Determines the length variations, and MNVs that cross splice junctions (both strands, acceptor and donor site). Lays foundation for a correct reporting of splice-junctions variants. The correct description of such variants requires reprting of coding region and intron alterations in a combined fashion (per HGSV guidelines).

8_x_POLYA:
Determines if polyA candidate sites at the 3' region are mutated. Currently, the AATAAA and ATTAAA and reverse complements thereof are interpreted as poly-adenylation sites. If needed, other consensus sites, not just a polyA, can be analyzed. The candidate mutations in polyA candidate sites should be evaluated using algorithms specifially designed for polyA-discovery. For example, DNA Functional Site Miner at http://dnafsminer.bic.nus.edu.sg. Same goes to any other domain that can be programmed in this script. Basically, any consenses DNA domaine witing known transcripts can be entered into this script.

9_x_CODING:
Splits a file harboring all the information from the upstream steps into the coding and non-coding region variant file. From this script on, only the coding region variants are analyzed.

10_x_T_MASTER:
Writes 24 chromosome specific scripts that fetch transcripts based on the coordinates provided by UCSC Genome Browser. Transcripts delived by these scripts will be spliced and translated in downstream steps as both Ref and Alt transcripts.

11_x_PRE_TRANS
Mutates all refence transctips as needed, and puts the mutated transcript into a 'Comparison' folder in each chromosome folder.

12_x_REF_TRANS.pl
Splices and translates each reference transcripts in silico. Works in combination with z_TRANSLATE.pl (see below).

13_x_Cx_COORDINATE.pl
Adjusts the trascript coordinates of each length variation. The variants crossing splice junction undergo a coordinate adjustment differnt from variants that only affect exons. If not done appropritely, the downstream variant translation and aligment with reference is inaccurate.

14_x_ALT_TRANS.pl
Translates the variant transcripts and alignes with the appropriate reference transcript from 12_x_REF_TRANS and z_TRANSLATE (see below). Works in combination with z_VARIANT_TRANSLATE.pl (see below).

15_x_PIPE_MASTER.pl
Runs the whole pipeline

z_TRANSLATE:
Code to carry out in silico splicing, and translation of reference transcripts

z_VARIANT_TRANSLATE.pl
Code to carry out in silico splicing, and translation of alt transcripts. Does Ref vs Alt alignments, both CDS and proteins, and prints a combined output.

CLEAN_UP.pl
Cleans the folders used in the pipeline. Does not touch the input and Archive folder.








