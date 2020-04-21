#! /usr/bin/env python

import os, sys, math, time
from optparse import OptionParser
from collections import defaultdict
import gzip

#Comai Lab, Ucdavis Genome Center
#Meric Lieberman, Isabelle Henry, 2019
# This work is the property of UC Davis Genome Center - Comai Lab

# Use at your own risk. 
# We cannot provide support.
# All information obtained/inferred with this script is without any 
# implied warranty of fitness for any purpose or use whatsoever. 
#------------------------------------------------------------------------------

#This script outputs a read coverage by bin across a reference sequence, using a directory of samtools aligned .sam files as input. 
#It can also output a measure of relative coverage compared to a control dataset. There can be two types of control data: either a 
#control file is indicated or the mean of all files in the directory is calculated and used as the control set. In both cases, the 
#values for relative percentage per bin were calculated by dividing the percentage of reads mapping to that bin for the sample at 
#hand by the mean percentage of reads mapping to that bin for the control set. Finally, all values are multiplied by the ploidy 
#parameter (default 2) such that values for bins present in X copies would oscillate around X.
#
#This script also outputs a second small file containing the number of read processed from each sam file.
#
#Usage: [...] denotes optional parameters, if not indicated, default parameters are used.
#bin-by-sam.py -o output-bin-file.txt -s size-of-bins [-c control .sam file] [-u] [-m number of max snps, default is 5] [-b] [-r] [-p ploidy for relative percent calculation] [-C]
#
#For help
#bin-by-sam.py -h
#
#Input:
#Run in a directory with the input .sam files. If you want to use one of the files as control for the relative coverage, specify the file with the -c option.
#
#Parameters
#
#Required:
#-o, output file name
#-s, bin size (bps)
#
#Optional, see below
#
#Output:
#One file with a line per bin of each reference sequence and a column for each input .sam library, as well as the relative coverage per input .sam library.

usage = "\n\n%prog"
parser = OptionParser(usage=usage)
parser.add_option("-c", "--controlfile", dest="f", default="NA", help="Input zipped sam file if wanting a specific library to be a control for normalization.")
parser.add_option("-o", "--out", dest="o", default = "NA", help="Output bin file.")
parser.add_option("-q", "--minqual", dest="minqual", type = "int", default=0, help="Min mapping quality")
parser.add_option("-s", "--binsize", dest="binsize", type = "int", default=1000000, help="Bin size")
parser.add_option("-b", "--breaks", dest="breaks", action="store_true", default = False, help="Insert breaks")
parser.add_option("-r", "--removefile", dest="r", default=False, help="A sam header file of reference sequences to ignore")
parser.add_option("-p", "--ploidy", dest="ploidy", type = "int", default=2, help="Ploidy multiplier, default is 2 for diploid")
parser.add_option("-C", "--covmode", dest="covmode", action="store_true", default = False, help="Only output coverage columns, not relative percent")
parser.add_option("-m", "--mode", dest="mode", type = "str", default = "-", help="Modes (S or TP are most common): \
    \t\t\t\t\tS - use this mode if reads are mapped single ended. \
    \t\t\t\t\tPS - reads are mapped paired, but count as if single ended. \
    \t\t\t\t\t\t\tTP - eads are mapped paired, this uses the default \"Correct\" PE mapping flags (For paired end use this unless understand the other options) \
    \t\t\t\t\t\t\tTPI - same as TP, but allow odd calculated inserts up to 2kb \
    \t\t\t\t\t\t\t\tTPM - same as TP, but allow same mapping direction for pairs -,- and +,+. \
    \t\t\t\t\t\t\t\tTPA - same as TP, but allow both odd inserts up to 2kb and same mapping direction for pairs -,- and +,+.")
parser.add_option("-l", "--listfile", dest="filelist", default="NA", help="A list of sam or sam.gz files to run on, this would be instead of runnign all sam files in current directory")
parser.add_option("-B", "--bamfile", dest="bams", action="store_true", default = False, help="Use unsorted .bam files instead of .sam or .sam.gz files. DO NOT USE SORTED BAMS!!!!")


(opt, args) = parser.parse_args()



outmode = False



if opt.mode not in ["S", "PS", "TP", "TPI", "TPM", "TPA"]:
   parser.error("Please specify a run mode with the -m paramter. S, PS, TP, TPI, TPM, or TPA. Use the -h option for paramter description")

if opt.o == "NA":
   parser.error("Please specify an output file using the -o parameter. Use the -h option for paramter description")

#parser.add_option("--unique", "-u", dest="unique", action="store_true", default = False, help="U/all only U (Off)")

if opt.f != "NA" and opt.covmode == True:
   parser.error("Cannot specify a contol then supress contol relative coverage percent columns")

#takes in a sam header file of chroms/genes to ignore, must be specified by command line
remlist = []
remcount = {}
remsize = {}
if opt.r != False:
   rem = open(opt.r)
   while 1:
      x = rem.readline()
      if x == '':
         break
      if x[0] != '@':
         break
      x[3333]
      if x[:3] == "@SQ":
         temp = (x[:-1].replace('SN:','').replace('LN:','').split('\t')[1:])
         key2 = temp[0]
         remlist.append(key2)
         remcount[key2] = 0
         remsize[key2] = int(temp[1])
   rem.close()

if opt.filelist != "NA":
   #read in list of sam files, must be in current directory and end in "_aln.sam"
   f = open(opt.filelist)
   todo = []
   for l in f:
      todo.append(l.replace('\n',''))
   f.close()
else:
   #read in list of sam files, must be in current directory and end in "_aln.sam"
   if opt.bams == False:
      li = os.listdir(os.getcwd())
      todo = filter(lambda x: ".sam" in x or ".sam.gz" in x, li)
   else:
      li = os.listdir(os.getcwd())
      todo = filter(lambda x: ".bam" in x, li)      



todo.sort()





#read sam header of chrom/genes to use
data = defaultdict(lambda : defaultdict(lambda: defaultdict(lambda: defaultdict(int))))
all = []
sizes = []
lookup = {}

if opt.bams == False:
   if todo[0].endswith('.gz'):
      f = gzip.open(todo[0], 'rb')
   else:
      f = open(todo[0])   
   while 1:
      x = f.readline()
      if x[0] != '@':
         break
      if "\tLN:0\n" in x:
         continue
      if x[:3] == "@SQ":
         temp = (x[:-1].replace('SN:','').replace('LN:','').split('\t')[1:])
         key2 = temp[0]
         if key2 in ["ChrUn", "ChrSy"]:
            continue
         if key2 not in remlist:
            all.append(key2)
            sizes.append(int(temp[1]))
            lookup[key2] = int(temp[1])
else:
   import pysam, gc
   f = pysam.AlignmentFile(todo[0], "rb")
   head = f.header
   head = head.to_dict()
   all2 = head['SQ']
   for item in all2:
      if type(item) == str:
         break
      tsize = int(item['LN'])
      ref = item['SN']
      if ref in ["ChrUn", "ChrSy"]:
         continue
      if ref not in remlist:
         all.append(ref)
         sizes.append(tsize)
         lookup[ref] = tsize
   f.close()
   del(all2)
   del(head)
   gc.collect()


#for data export purpose; add blank lines based on size of largest reference
numblanks = max(sizes)/opt.binsize/10

fseen = []

f.close()
globalcount = {}
count = 0
liblist = []
if opt.mode == "S":
   #per sam file, count reads
   for file in todo:
      if file.endswith(".bam") == False:
         if file.endswith(".gz"):
            f = gzip.open(file, 'rb')
         else:
            f = open(file)
      else:
         f = pysam.AlignmentFile(file, "rb")
      print file
      libname = file.split('.')[0].replace('_aln','')
      liblist.append(libname)
      globalcount[libname] = 0
      count = 0
      #print time.asctime()
      for x in f:
         if opt.bams == True:
            x = x.tostring()
         count +=1
         if count % 1000000 == 8:
            print count
   
         if x[0] == '@':
            continue       
         
         l = x.replace('\n','').split('\t')   
         if l[2] == '*':
            continue
         if int(l[1]) > 80 and int(l[1]) < 2000:
            parser.error("Please specify a correct run mode - S, PS, TP, TPI, TPM, or TPA. Use the -h option for paramter description")
            break
         if l[1] in ['0', '16']:
            pos = l[3]   
         else:
            continue

         if int(l[4]) < opt.minqual:
            continue

         key1 = l[2]
         key2 = l[2]
         if key2 in remlist:
            remcount[key2] += 1
            continue
        
         s1 =  int(l[3])
         e1 = s1 + len(l[9])

         mid = (e1 - s1 +1) / 2 + s1
         key3 = int(mid) / opt.binsize
         data[key1][key2][key3][libname] += 1
         globalcount[libname] +=1
      f.close()
else:
   #per sam file, count reads
   
   for file in todo:
      if file.endswith(".bam") == False:
         if file.endswith(".gz"):
            f = gzip.open(file, 'rb')
         else:
            f = open(file)
      else:
         f = pysam.AlignmentFile(file, "rb")      
      print file
      libname = file.split('.')[0].replace('_aln','')
      liblist.append(libname)
      globalcount[libname] = 0
      count = 0
      #print time.asctime()
      while 1:
         if opt.bams == False:
            x1 = f.readline()
         else:
            try:
               xt = f.next()
            except StopIteration:
               break
            x1 = xt.tostring() 
         count +=1
         if count % 1000000 == 8:
            print count
   

         if x1 == '':
            break
         if x1[0] == '@':
            continue 
         l1 = x1.replace('\n','').split('\t')
         while int(l1[1]) > 200:
            if opt.bams == False:
               x1 = f.readline()
            else:
               xt = f.next()
               x1 = xt.tostring() 
            l1 = x1.replace('\n','').split('\t')
               
               
         if opt.bams == False:
            x2 = f.readline()
         else:
            xt2 = f.next()
            x2 = xt2.tostring() 
         l2 = x2.replace('\n','').split('\t')         
         while int(l2[1]) > 200:
            if opt.bams == False:
               x2 = f.readline()
            else:
               xt2 = f.next()
               x2 = xt2.tostring() 
            l2 = x2.replace('\n','').split('\t')
   
         #look for slipped pair
         if int(l1[1]) > 200 or int(l2[1]) > 200:
            x[123312321132123132]
         #look for slipped pair
         if l1[0] != l2[0]:
            x[372678231]
         
         flags = [int(l1[1]), int(l2[1])]
         flags.sort()
         if flags[0] < 60:
            parser.error("Please specify a correct run mode - S, PS, TP, TPI, TPM, or TPA. Use the -h option for paramter description")

         
         if flags in [[99, 147], [83, 163]]:
            #correctly normal mapped
            pass
         elif flags in [[81, 161],[97, 145]] and abs(int(l1[8])) < 2001 and opt.mode in ["TPI", "TPA", "PS"]:
            #weird insert, but mapped in correct orientations, and in an allowing mode
            pass
         elif flags in [[67, 131],[115, 179]] and opt.mode in ["TPM", "TPA","PS"]:
            # mappen in wrong orientation -,- or +,+, but good otherwise
            pass
         elif flags in [[65, 129],[113, 177]] and abs(int(l1[8])) < 2001 and opt.mode in ["TPA", "PS"]:
            #mapped in wrong orientation and weird insert
            pass
         elif flags == [77,141]:
            continue
         elif opt.mode == "PS":
            #print l1, '\n', l2, '\n\n'
            pass  
         else:
            continue
            
         if flags not in fseen:
            fseen.append(flags)
            print flags
            

         
         
         s1 =  int(l1[3])
         e1 = s1 + len(l1[9])
         s2 = int(l2[3])
         e2 = s2 + len(l2[9])
         zone = [min(s1, e1, s2, e2), max(s1, e1, s2, e2)]
         midpt = ((zone[1] - zone[0]+1) / 2) + zone[0]
         
         if opt.mode != "PS":
            if int(l1[4]) < opt.minqual or int(l2[4]) < opt.minqual:
               continue      
            key1 = l1[2]
            key2 = l1[2]
            if key2 not in remlist:
               key3 = int(midpt) / opt.binsize         
               data[key1][key2][key3][libname] += 1
               globalcount[libname] +=1
            else:
               remcount[key2] += 1
         else:
            if l1[5] != '*':
               if int(l1[4]) > opt.minqual:            
                  key1 = l1[2]
                  key2 = l1[2]
                  if key2 in remlist:
                     remcount[key2] += 1
                  else:
                     key3 = (s1+(e1-s1+1)/2) / opt.binsize         
                     data[key1][key2][key3][libname] += 1
                     globalcount[libname] +=1
            
            if l2[5] != '*':
               if int(l2[4]) > opt.minqual:  
                  key1 = l2[2]
                  key2 = l2[2]
                  if key2 in remlist:
                     remcount[key2] += 1
                  else:
                     key3 = (s2+(e2-s2+1)/2) / opt.binsize         
                     data[key1][key2][key3][libname] += 1
                     globalcount[libname] +=1
      f.close()

#use control lib if specified, otherwise non applicable
control = opt.f
control = control.split('.')[0].replace('_aln','')

#create header for output file
header = ['Chrom', 'Strt', 'End']
header+= liblist
if opt.covmode == False:
   header+= map(lambda x: x+"/"+control, liblist)
o = open(opt.o, 'w')
o.write('\t'.join(header)+'\n')

#traverse each chromosome by bin
#then in each bin, first tabulate relative %
#then apply to create additional columns for relative%/control(or all)%
#when have one set for cov and other set for %, output line and repeat
#will output blank lines as/if speciifed
# res = defaultdict(list)
# for chrom in all:
#    part = chrom[:]
#    bins = data[part][chrom].keys()
#    bins.sort()
#    for modbin in bins:
#       libdata = data[part][chrom][modbin]             
#       line = [chrom, modbin*opt.binsize+1, (modbin+1)*opt.binsize]
#       perst = {}
#       for x in liblist:
#          try:
#             temper = libdata[x]/float(globalcount[x])
#          except:
#             temper = 0.0
#             
#          res[x].append([temper, chrom, modbin])
# 
# badbins = []
# if outmode == True:
#    for lib in res:
#       hold = res[lib]
#       alldat = map(lambda z: z[0], hold)
#       mn = np.mean(alldat)
#       sd = np.std(alldat)
#       #iqr = np.percentile(dist, 75) - np.percentile(dist, 25)
#       bad = filter(lambda x: x[0]<(mn-outmult*sd) or x[0] > (mn+outmult*sd), hold)
#       if len(bad) > 0:
#          for item in bad:
#             ct = item[1]
#             mt = item[2] 
#             if mt in data[ct][ct]:
#                del data[ct][ct][mt]
#                badbins.append(ct+'_'+str(mt*opt.binsize+1)+'_'+str((mt+1)*opt.binsize))



for chrom in all:
   part = chrom[:]
   bins = data[part][chrom].keys()
   bins.sort()
   for modbin in bins:
      libdata = data[part][chrom][modbin]             
      line = [chrom, modbin*opt.binsize+1, (modbin+1)*opt.binsize]
      if modbin == bins[-1]:
         line = [chrom, modbin*opt.binsize+1, lookup[chrom]]
      pers = {}
      for x in liblist:
         if libdata[x] != 0:
            line.append(libdata[x])
            pers[x] = libdata[x]/float(globalcount[x])
         else:
            line.append(0)
            pers[x] = 0.0
      sums = sum(pers.values())/float(len(pers.values()))
      if opt.covmode == False:                  
         for x in liblist:
            if control == 'NA':
               try:
                  line.append(round(pers[x]/sums*opt.ploidy, 3))
               except ZeroDivisionError:
                  line.append(0.0)
            else:
               try:
                  line.append(round(pers[x]/pers[control]*opt.ploidy, 3))
               except ZeroDivisionError:
                  line.append('.')
      fline = map(lambda x: str(x), line)
      o.write('\t'.join(fline)+'\n')
   if opt.breaks == True:
      o.write(''.join(map(lambda x: '\n', range(numblanks))))

o.close()

#file to hold library statistics
#and output by each lib
o2 = open('readcounts-'+opt.o, 'w')

sechead = ['Lib', 'Reads', 'Reads/MB']
o2.write('\t'.join(sechead)+'\n')
tot = sum(sizes)/1000000.0
for x in liblist:
   o2.write('\t'.join([x, str(globalcount[x]), str(round(globalcount[x]/tot, 2))])+'\n')
if opt.r != False:
   o2.write('\n\nRemoved Reference Counts:\n')
   o2.write('Reference\tReads\n')
   for x in remlist:
      o2.write(x+'\t'+str(remcount[x])+'\n')

o2.close()           
         
         
      
   
   
   
   
   
   
      
      
