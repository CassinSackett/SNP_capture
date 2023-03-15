# Population genomics pipeline from SNP data with a reference genome

###### tags: `SNPs`, `Capture`, `in-solution hybridization`, `GATK`


> This pipeline is designed for modern DNA of non-model organisms but does use a reference genome for alignment :bird: 

- Table of Contents
[ToC]

## Part I: Read filtering, assembly, quality control

### Step 1: Check the quality of samples and trim/filter accordingly

#### 1a. Check the md5 for each file to make sure it downloaded/uploaded correctly.
```
md5sum sample1_R1.fastq.gz
```
or, to run in a loop (best done within a batch script for multiple large files; see step 1b for detailed structure of a batch script):
```
#!/bin/bash
#SBATCH ...
#SBATCH ...

for i in *.fastq.gz; do md5sum $i; done
```
On some computers, the command will be ```md5``` instead of ```md5sum```.

Copy the md5 codes from your terminal output and compare them to the files sent by the sequencing company. I use the lazy way and copy all the md5 info into one text file and open it in BBEdit. Then when you put the cursor within one of the codes, it will underline the whole code if that sequence appears again in the file.

#### 1b. Run FastQC on compressed files - this may take a long time to uncompress the files
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

#### 1c. Copy html files to your local drive to open them in a browser.  
In a shell for your local machine, navigate to where you want the files and type 
```javascript    
scp yourname@qbc.loni.org:/work/yourname/*html ./
```

#### 1d. Based on fastqc results, trim files to a length that contains a large proportion of high-quality reads (usually Q20+). 
I run Trimmomatic in paired-end mode, but it will keep unpaired reads if you want to use them later.
```javascript
#SBATCH lines
#SBATCH -e trimmomatic_clipNextera_031023.err

module load jdk/1.8.0_262/intel-19.0.5

for i in /work/yourname/*1.fq; do java -jar /project/sackettl/Trimmomatic-0.39/trimmomatic-0.39.jar $i ${i%1.fq}2.fq -baseout ${i%1.fq} ILLUMINACLIP:/project/sackettl/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:60; done
```

#### 1e. Re-run fastqc on the trimmed samples and compare the number of reads before and after quality filtering. 
If you have a ton of samples, you can just look at a few to get an idea of how much changed. After this, use a batch script to compress original fq files to save space (``` for i in *1.fastq; do gzip $i; done```).

### Step 2: Align reads to your reference genome
#### 2a. If your reference genome has not been indexed (i.e., if the only file in the directory is the .fasta file), you should index it now.  
First, create a directory for all the reference files, move the .fasta file there, and cd to that directory. Then run
```
#!/bin/bash
#SBATCH ...
#SBATCH ...

bwa index -a bwtsw genome.fasta 

samtools faidx genome.fasta
```
where ```bwa index``` will create *.bwt, *.pac, *.ann, *.amb and *.sa files. 

and ```samtools faidx```  will create a *.fai file. All of these are necessary for the next step. Later you’ll need GATK to create a *.dict file for the reference, or you can do it now with 
```
home/gatk-4.1.2.0/gatk --java-options "-Xmx2G" CreateSequenceDictionary -R /sackettl/ref.fasta
```

#### 2b. Now you are ready to align your sequences!  :elephant:
You will do so running BWA mem on the filtered fastq files (bwa mem is an algorithm for long reads (> 100bp) and split alignment; it is recommended for high-quality input sequences because it is more accurate than some other aligners). For paired end reads, reads 1 and 2 must be in the same order. If you used the ```fastq_trimmer.py``` script to trim, your script will look like this: 
```
#!/bin/bash
#SBATCH ...

for i in *R1_trim.fastq; do bwa mem -v 3 -M -P -a -t 10 /ref/ref.fasta $i $
{i/R1trim/R2trim} > ${i%.fastq}.sam 2> ${i%.fastq}.mem.log; done
```
or if you used trimmomatic, it will look more like this:
```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

module load bwa/0.7.17/intel-19.0.5

for i in /pdog/trimmed/*1P; do bwa mem -v 3 -M -P -a -t 44 /sackettl/pdogkrakengenome/pilon_pdog_kraken_nothuman_1308.fasta $i ${i%1P}2P > ${i%1P}PE.sam 2> ${i%1P}PE.mem.log; done
```

which produces a single paired-end output file for each pair of inputs. This is mapping only the reads that contained pairs after Trimmomatic's filtering (you can use the others if you want them; they have filenames ending in U for Unpaired.
Parameter flags should be listed in this order (v first), where 
* –M: mark shorter split hits as secondary (important for Picard compatibility/ functionality with MarkDuplicates)
* –t: # threads 
* –P: In paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair 
* –a: Output all alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.

### Step 3: Convert .sam files to bam files for downstream analyses
If you mapped paired reads and did not quality filter bases/reads previously (i.e., if you went straight from Trimmomatic to bwa), do so now by adding the q flag (before the other flags) and a mapping quality threshold, e.g., ```-q 30```
```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

for i in *.sam; do samtools view -q 20 -bt /pdogkrakengenome/pilon_pdog_kraken_nothuman_1308.fasta -o ${i%sam}bam $i; done 
```

### Step 4: Process bam files and prepare them for downstream processing

#### 4a. Add read group IDs and sort the bam files 


```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

cd $PBS_O_WORKDIR
mkdir -p tmp

source activate picard

for i in /pdog/*PE.bam; 
do 
        java -Dpicard.useLegacyParser=false -Xmx2g -jar $CONDA_PREFIX/share/picard-2.26.11-0/picard.jar \
        AddOrReplaceReadGroups INPUT=$i \ 
        OUTPUT=${i%_PE.bam}.tag.bam MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$PWD/tmp \
        SO=coordinate RGID=${i%_PE.bam} RGLB=1 RGPL=illumina RGPU=1 RGSM=${i%_PE.bam}; 
done 

```

where
* RGID is Read Group ID. Currently this is set to the individual, but you could do modern vs ancient or high/low/mid (e.g., in sample names, every second underscore is the RGID). This modifies the bam headers so that they are “tagged” with the read group ID (necessary for downstream analyses).  
* MAX_RECORDS_IN_RAM is necessary for large files (e.g., over 60 GB)
* TMP_DIR is necessary if there are limits to the number of files allowed in a directory
* SO = sort order; sorts by where the reads align in the reference genome

Many enormous temp files are created in this step, so if you install yourself on a HPC, you need to add a tmp directory in a location with unlimited storage so they won’t be saved to the compute node. Before my command, I added ```cd $PBS_O_WORKDIR \ mkdir -p tmp``` and then in my command before the SO part I added ```TMP_DIR=$PWD/tmp``` to solve this problem.

Files can have only one tag, so if you want to tag by some group membership, you may have to do this twice – once for lane/run number before merging, and once for your tag of interest after the replicates are merged.

If you have trouble running this step, it may be because of the way picard tools is set up (it is not C-compiled, so you cannot add it to your path -- you have to set an environment variable). You can set up an alias called runpicard and run it like this:
```
i in *.bam; do runpicard AddOrReplaceReadGroups ...
```


#### 4b. If samples were sequenced on multiple runs or lanes, merge reads now. 
 
All of the .bam files need to be in the same folder for this step.
```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

cd $PBS_O_WORKDIR
mkdir -p tmp

source activate picard

for i in /pdog/*L001.tag.bam; 
do 
        java -Dpicard.useLegacyParser=false -Xmx2g -jar $CONDA_PREFIX/share/picard-2.26.11-0/picard.jar \
        MergeSamFiles INPUT=$i INPUT=${i%L001_tag.bam}L002_tag.bam INPUT=${i%L001_tag.bam}L007_tag.bam INPUT=${i%L001_tag.bam}L008_tag.bam \
        OUTPUT=${i%L001_tag.bam}.merged.bam MAX_RECORDS_IN_RAM=1000000 TMP_DIR=$PWD/tmp;
done 

```

#### 4c. Remove PCR duplicates
This step removes PCR duplicates. You should look at the metrics file for statistics on how many were duplicates, etc. 

```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

cd $PBS_O_WORKDIR
mkdir -p tmp

source activate picard

for i in /project/sackettl/MolEvol/pdog/*merged.bam;  
do
        java -Dpicard.useLegacyParser=false -Xmx2g -jar $CONDA_PREFIX/share/picard-2.26.11-0/picard.jar \
        MarkDuplicates INPUT=$i OUTPUT=${i%merged.bam}rmdup.bam MAX_RECORDS_IN_RAM=1000000 \
        METRICS_FILE=${i%merged.bam}rmdup.metrics ASSUME_SORTED=true;
done

```

#### 4d. Index sorted, duplicate-filtered bam files

All bam files (.tag.bam, etc.) need to be in the same folder for this step.

```
#!/bin/bash
#SBATCH -n 48
#SBATCH ...

for i in *.rmdup.bam; do samtools index $i; done
```

Now you are ready to genotype your samples! :boom: :boom:


## Part II: Genotyping




## Part III: Population Genomic Analyses



