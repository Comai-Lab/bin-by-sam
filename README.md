# bin-by-sam

Comai Lab, Ucdavis Genome Center
Meric Lieberman, Isabelle Henry, 2019
This work is the property of UC Davis Genome Center - Comai Lab

Use at your own risk. 
We cannot provide support.
All information obtained/inferred with this script is without any 
implied warranty of fitness for any purpose or use whatsoever. 

------------------------------------------------------------------------------

This script outputs a read coverage by bin across a reference sequence, using a directory of samtools aligned .sam files as input. 
It can also output a measure of relative coverage compared to a control dataset. There can be two types of control data: either a control file is indicated or the mean of all files in the directory is calculated and used as the control set. In both cases, the values for relative percentage per bin were calculated by dividing the percentage of reads mapping to that bin for the sample at hand by the mean percentage of reads mapping to that bin for the control set. Finally, all values are multiplied by the ploidy parameter (default 2) such that values for bins present in X copies would oscillate around X.

This script also outputs a second small file containing the number of read processed from each sam file.

Usage: [...] denotes optional parameters, if not indicated, default parameters are used.
bin-by-sam.py -o output-bin-file.txt -s size-of-bins [-c control .sam file] [-u] [-m number of max snps, default is 5] [-b] [-r] [-p ploidy for relative percent calculation] [-C]

For help
bin-by-sam_v7.py -h

Input:
Run in a directory with the input .sam (or .bam or sam.gz) files. If you want to use one of the files as control for the relative coverage, specify the file with the -c option.

Parameters

Required:
-o, output file name
-s, bin size (bps)
-m, mode - this is the read type. (S or TP are most common) The options are:
   S - use this mode if reads are mapped single ended.
   PS - reads are mapped paired, but count as if single ended.
   TP - reads are mapped paired, this uses the default "Correct" PE mapping flags (For paired end use this unless understand the other options)
   TPI - same as TP, but allow odd calculated inserts up to 2kb
   TPM - same as TP, but allow same mapping direction for pairs -,- and +,+.
   TPA - same as TP, but allow both odd inserts up to 2kb and same mapping direction for pairs -,- and +,+.

Optional:
"-c" or "--controlfile", This is for an input sam specific library to be a control for normalization.
"-q" or "--minqual", This is for the minimum mapping quality
"-b" or "--breaks", This puts several lines of sapcing between references.
"-r" or "--removefile", A sam file header of sequences not to use 
"-p" or "--ploidy", Ploidy multiplier, default is 2 for diploid
"-C" or "--covmode", This is to only output coverage columns, not the relative percent columns as well
"-B" or "--bamfile", This uses unsorted .bam files instead of .sam or .sam.gz files. DO NOT USE SORTED BAMS!

Output:
One file with a line per bin of each reference sequence and a column for each input .sam library, as well as the relative coverage per input .sam library.
