#!/usr/bin/python

import sys
import os
import numpy
import csv

fstfile = open('amakihiR1R2_DP9mendelrm25pctind80pctloc65pctind90pctloc_PuuLaau-Pahoa.weir.fst', 'r')
outfile = open('PuuLaau-Pahoa_highest-fst.txt', 'w')
fstfilereader = csv.reader(fstfile, delimiter='\t')
highfst = csv.writer(outfile, delimiter='\t'

for row in fstfilereader:
	if row['WEIR_AND_COCKERHAM_FST'] > 0.9:  #i want this to say if 3rd element is greater than x, print line
		highfst.writerow(row)







# Write results to files
#print "writing file"
#numpy.savetxt('PuuLaau-Pahoa_highest-fst.txt', outfile, fmt='%s', delimiter="\t")


