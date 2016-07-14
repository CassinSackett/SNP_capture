#!/usr/bin/python

# I will need to find lines beginning with x, split the line into components, and extract the ith component.
# Then I will have to store that component in a file, along with some information about which file it came from
# OR other text in that file.

# The input file is space delimited and does not have a format common to all lines
# mean value of ln likelihood can be found on line 6158
# Estimated Ln Prob of Data can be found on line 6157
# K can be found in line 1 (qsub command), 50, 6143, 6150 ("3 populations assumed"), 6179=last line ("job xyz done")
# Unfortunately, almost everything over K=4 terminated prematurely :(

#import statements
import sys


# Find the line with the name of the job script or output file (with K information)
def GetKvalue(line):

	#split on whitespace
	Kdata = line.split()
	
    # check the length of the line (we obviously need >=4 elements if we're going to operate on the 4th element)
	if len(Kdata) >= 2:

        # Make sure the desired element is an integer  
		try:
			Kval = int(Kdata[0])
		except ValueError:	
			pass

    # return the value of K
	return Kdata[0]


# split the line with the likelihood value on the '=' sign and keep what's after the =
def GetLL(line):

    # split on tabs
    likelihooddata = line.split('=')

    # check the length of the line (we obviously need 2 elements if we're going to operate on the 2nd element)
    if len(likelihooddata) >= 2:

        # Make sure the second element is either an integer or a float 
        try:
            lnL = int(likelihooddata[1])
        except ValueError:
            pass
        try:
            lnL = float(likelihooddata[1])
        except ValueError:
            pass

    # print the value of the log likelihood
    return likelihooddata[1]

# tell python "the program starts here". 
# The program starts executing near the bottom because you have to have all your functions 'seen' by the interpreter before you can call them.
# All functions or classes have to be above this point.
if __name__ == "__main__":

#    filename = "structure_Hawaiiall2famgrp_K3_013116.log"
	filename = sys.argv[1]
	outfile = open(filename + 'LL', 'w') #it would be better to extract info from the filename and rename accordingly...

	with open(filename,"r") as logFile:
		for line in logFile:
			if line.endswith('populations assumed\n'):
				line = GetKvalue(line)
				outfile.write(line[0:2]) #How can I insert a tab between this and the next values? 
			elif line.startswith('Mean value of ln likelihood'):
				line = GetLL(line)
				outfile.write(line[1:66]) #I'm too lazy to figure out how to say 1 until the end of the line

# I need some way to join these 2 values with a tab instead of printing them with no delimitation


