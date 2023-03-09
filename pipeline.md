# Population genomics pipeline from SNP data with a reference genome

###### tags: `SNPs`, `Capture`, `in-solution hybridization`, `GATK`


> This pipeline is designed for modern DNA of non-model organisms but does use a reference genome for alignment :bird: 

- Table of Contents
[ToC]

## Part 1: Read filtering, assembly, quality control

### Step 1: Uncompress files, check the quality of samples, and trim accordingly
- 1a. Check the md5 for each file to make sure it downloaded/uploaded correctly.
```
md5 sample1_R1.fastq.gz
```
or, to run in a loop (best done within a batch script for large files; see step 1b for structure of a batch script):
```
for i in *.fastq.gz; do md5 $i; done
```


- 1b. Run FastQC on compressed files - this may take a long time to uncompress the files
```javascript
#!/bin/bash
#SBATCH -q workq
#SBATCH -N 1
#SBATCH -n 18
#SBATCH -t 6:00:00
#SBATCH -A loni_molec_evo
#SBATCH -o fastqc_030823.out
#SBATCH -e fastqc_030823.err

/project/sackettl/FastQC/fastqc *fastq.gz --threads 16 --outdir=/work/yourname/
```

- 1b. Copy html files to your local drive to open them in a browser.  In a shell for your local machine, navigate to where you want the files and type 
```javascript    
scp yourname@qbc.loni.org:/work/yourname/*html ./
```

- 1c. Based on fastqc results, trim files to a length that contains a large proportion of high-quality reads (usually Q20+). I run Trimmomatic in paired-end mode, but it will keep unpaired reads if you want to use them later.
```javascript
#SBATCH lines
#SBATCH -e trimmomatic_clipNextera_031023.err

module load jdk/1.8.0_262/intel-19.0.5

for i in /work/yourname/*1.fq; do java -jar /project/sackettl/Trimmomatic-0.39/trimmomatic-0.39.jar $i ${i%1.fq}2.fq -baseout ${i%1.fq} ILLUMINACLIP:/project/sackettl/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:60; done
```

- 1d. Re-run fastqc on the trimmed samples and compare the number of reads before and after quality filtering. In the future, there will be a Python script to extract these numbers automatically, but you can just look at a few. After this, use a batch script to compress original fq files to save space (``` for i in *1.fastq; do gzip $i; done```).

### Step 2: Align reads to your reference genome
- 2a. If your reference genome has not been indexed (i.e., if the only file in the directory is the .fasta file), you should index it now.  First, create a directory for all the reference files, move the .fasta file there, and cd to that directory. Then run
```
bwa index -a bwtsw genome.fasta 
```
which will create *.bwt, *.pac, *.ann, *.amb and *.sa files. 

Next, run
```
samtools faidx genome.fasta
```
which will create a *.fai file. All of these are necessary for the next step. Later you’ll need GATK to create a *.dict file for the reference, or you can do it now with 
```
home/gatk-4.1.2.0/gatk --java-options "-Xmx2G" CreateSequenceDictionary -R /sackettl/ref.fasta
```

- 2b. Now you are ready to align your sequences!  :elephant:
You will do so running BWA mem on the filtered fastq files (bwa mem is an algorithm for long reads (> 100bp) and split alignment; it is recommended for high-quality input sequences because it is more accurate than some other aligners). For paired end reads, reads 1 and 2 must be in the same order.
```
for i in *R1_trim.fastq; do bwa mem -v 3 -M -P -a -t 10 /ref/ref.fasta $i $
{i/R1trim/R2trim} > ${i%.fastq}.sam 2> ${i%.fastq}.mem.log; done
```
which produces a single paired-end output file for each pair of inputs. Parameter flags should be listed in this order (v first), where 
* –M: mark shorter split hits as secondary (important for Picard compatibility/ functionality with MarkDuplicates)
* –t: # threads 
* –P: In paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair 
* –a: Output all alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.




## Part 2: Genotyping




## Part 3: Population Genomic Analyses



