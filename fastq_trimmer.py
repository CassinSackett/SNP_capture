#!/usr/bin/python

# Program goes through a FastQ file and trims from 75 onward (removes last 25 bp)

# FastQ format:
# 	Line 1: Sequence ID
#	Line 2: Sequence
#	Line 3: '+' (and repeat of sequence ID)
#	Line 4: Quality scores

# import statements
import string
import sys

# define functions
def trimmer(line):
	if len(line) == 101:
		x = line[:-26]		# remove last 25 letters of line 
		return x + "\n"
	else:
		return line
	
if len(sys.argv)!= 2:
	print "USAGE: "+sys.argv[0]+" <FILENAME>"
	sys.exit()
	
# first, get input file from the user
#filename = raw_input('Enter the FastQ file: ')
filename = sys.argv[1]


# make the file available to use
file   = open(filename, 'r')
file2  = open(filename + 'NEW', 'w')
	
for line in file:
	file2.write(trimmer(line))
	