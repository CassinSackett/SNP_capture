# The purpose of this script is to remove all bad SNPs from a Lositan output file.
# The initial output file (here snpfile) file has the format SNP_number '\t' Heterozygosity '\t' Fst '\t' Pvalue
# But it also has Fst values of -100 (a sort of 'missing data' code) which can lead to spurious results in downstream calculations.
# We need to get rid of these values first.

snpfile = open('locuslist_nohead.txt')  
no_neg100sfile = open('locuslist_nohead_noneg100values.txt', 'w')  
outfile = open('significant_fstloci.txt', 'w')


for line in snpfile:
	if '-100' in line:
		del line
	elif 'NA' in line:
		del line
	else: 
		no_neg100sfile.write(line)

# We may also want to look only at loci that have significant differentiation or some other criteria.
no_neg100sfile = open('locuslist_nohead_noneg100values.txt', 'r')

for line in no_neg100sfile:
	columns = line.split('\t')
	pval = float(columns[3])
	if pval < 0.05:
		outfile.write(line)


