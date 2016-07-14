#!/usr/bin/python

outfile = open('amakihi_65pctind90pctloc_negsconverted.weir.fst', 'w')

# The input file format (tab delimited) is: Chromosome	Position	FST
# a function to modify the 3rd column to be 0 if the value is a number and is less than 0
def ChangeNegativeFst(line):

    # split on tabs
    SNPdata = line.split('\t')

    # check the length of the line (we obviously need 3 elements if we're going to operate on the 3rd element)
    if len(SNPdata) >= 3:

        # Make sure the third element is either an integer or a float 
        try:
            FST = int(SNPdata[2])
        except ValueError:
            pass
        try:
            FST = float(SNPdata[2])
        except ValueError:
            pass
        # the question 'is it less than zero' is asked and if so, the value of the 3rd word is replaced with '0'.
        try:
            if FST <= 0:
                SNPdata[2] = "0\n"
        except NameError:
            pass

    # return the line with a tab separating each word
    return '\t'.join(SNPdata)

# tell python "the program starts here". 
# The program starts executing near the bottom because you have to have all your functions 'seen' by the interpreter before you can call them.
# All functions or classes have to be above this point.
if __name__ == "__main__":

    filename = "amakihi_65pctind90pctloc.weir.fst"

    with open(filename,"r") as theFile:
        for line in theFile:
            line = ChangeNegativeFst(line)
            outfile.write(line)
