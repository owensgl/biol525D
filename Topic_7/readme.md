---
title: "Topic 7: SNP calling with GATK"
permalink: /Topic_7/
topickey: 7
topictitle: "SNP Calling"
---

## Accompanying material
* [Slides](./Topic 7.pdf)

In this tutorial we're going to call SNPs with GATK. 

The first step is again to set up directories to put our incoming files.
```bash

cd ~
mkdir log
mkdir gvcf
mkdir db
mkdir vcf
```
We also have a few programs we're going to use. Since we will be calling them repeatedly, its helpful to save their full path to a variable. This will only last for the current session so if you log out you'll have to set them up again.
```bash
gatk=/mnt/bin/gatk-4.1.2.0/gatk
picard=/mnt/bin/picard.jar
```

There are 10 different samples and we're going to have to run multiple steps on each. To make this easier, we make a list of all the sample names.
```bash
ls bam | grep .sort.bam$ | sed s/.sort.bam//g > samplelist.txt
```
Lets break this down. 

**ls bam** <= List all the files in the _bam_ directory

**\| grep .sort.bam$** <= Only keep the file names ending in _.sort.bam_.

**\| sed s/.sort.bam//g** <= Replace the string _.sort.bam_ with "", effectively leaving only the sample name.

**> samplelist.txt** <= Save the result in a file named _samplelist.txt_


The first step is to make duplicate reads using picardtools. If you were using GBS data you wouldn't want to do this step.

```bash

while read name; do
  java -jar $picard MarkDuplicates \
  I=bam/$name.sort.bam O=bam/$name.sort.dedup.bam \
  M=log/$name.duplicateinfo.txt
  samtools index bam/$name.sort.dedup.bam
done < samplelist.txt

```

Now in the bam files duplicate reads are flagged. Take a look in the log directory, which sample has the highest number of duplicate reads?


To use GATK, we have to index our reference genome. An index is a way to allow rapid access to a very large file. For example, it can tell the program that the third chromosome starts at bit 100000, so when the program wants to access that chromosome it can jump directly there rather than scan the whole file. Some index files are human readable (like .fai files) while others are not.
```bash
java -jar $picard CreateSequenceDictionary R= ref/HanXRQr1.0-20151230.1mb.fa O= ref/HanXRQr1.0-20151230.1mb.dict

samtools faidx ref/HanXRQr1.0-20151230.1mb.fa
```
Take a look at the ref/HanXRQr1.0-20151230.1mb.fa.fai. How many chromosomes are there and how long is each? 



The next step is to use GATK to create a GVCF file for each sample. This file summarizes support for reference or alternate alleles at all positions in the genome. It's an intermediate file we need to use before we create our final VCF file.

This step can take a few minutes so lets first test it with a single sample to make sure it works.
```
for name in `cat ~/samplelist.txt | head -n 1`
do
$gatk --java-options "-Xmx15g" HaplotypeCaller \
   -R ref/HanXRQr1.0-20151230.1mb.fa \
   -I bam/$name.sort.dedup.bam \
   --native-pair-hmm-threads 3 \
   -ERC GVCF \
   -O gvcf/$name.sort.dedup.g.vcf
done
```
 Check your gvcf file to make sure it has a .idx index file. If the haplotypecaller crashes, it will produce a truncated gvcf file that will eventually crash the genotypegvcf step. Note that if you give genotypegvcf a truncated file without a idx file, it will produce an idx file itself, but it still won't work.

We would run the HaplotypeCaller on the rest of the samples, but that will take too much time, so once you're satisfied that your script works, you can copy the rest of the gvcf files (+ idx files) from /mnt/data/gvcf into ~/gvcf.



The next step is to import our gvcf files into a genomicsDB file. This is a compressed database representation of all the read data in our samples. It has two important features to remember:
1. Each time you call GenomicsDBImport, you create a database for a single interval. This means that you can parallelize it easier, for example by calling it once per chromosome.
2. The GenomicsDB file contains all the information of your GVCF files, but can't be added to, and can't be back transformed into a gvcf. That means if you get more samples, you can't just add them to your genomicdDB file, you have to go back to the gvcf files.


We need to create a map file to GATK where our gvcf files are and what sample is in each. Because we use a regular naming scheme for our samples, we can create that using a bash script.
This is what we're looking for:

sample1 \t gvcf/sample1.g.vcf.gz

sample2 \t gvcf/sample2.g.vcf.gz

sample3 \t gvcf/sample3.g.vcf.gz

```bash

for i in `ls gvcf | grep ".g.vcf" | grep -v ".idx" | sed s/.sort.dedup.g.vcf//g`
do
  echo -e "$i\tgvcf/$i.sort.dedup.g.vcf"
done > ~/biol525d.sample_map

```

Lets break down this loop to understand how its working 

**for i in `...`** <= This is going to take the output of the commands in the `...` and loop through it line by line, putting the each line into the variable $i. 

**ls gvcf** <= List all files in the gvcf directory.

**\| grep ".g.vcf"** <= Only keep the files including .g.vcf in their name.

**\| grep -v ".idx"** <= Remove any file including .idx, which are the index files

**\| sed s/.sort.dedup.g.vcf//g** <= Remove the suffix to the filename so that its only the sample name remaining.

**do** <= Starts the part of the script where you put the commands to be repeated.

**echo -e "$i\t$gvcf/$i.sort.dedup.g.vcf"** <= Print out the variable $i (which is the sample name) and then a second column with the full file name.

**done > ~/biol525d.sample_map** <= Take all the output and put it into a file name _biol525d.sample_map_.


Next we call GenomicsDBImport to actually create the database.
```bash
$gatk --java-options "-Xmx10g -Xms10g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path db/HanXRQChr01 \
       --batch-size 50 \
       -L HanXRQChr01 \
       --sample-name-map ~/biol525d.sample_map \
       --reader-threads 3
```

With the genomicsDB created, we're finally ready to actually call variants and output a vcf
```bash
$gatk --java-options "-Xmx10g" GenotypeGVCFs \
   -R ref/HanXRQr1.0-20151230.1mb.fa \
   -V gendb://db/HanXRQChr01 \
   -O vcf/HanXRQChr01.vcf.gz
```
Now we can actually look at the VCF file

```bash
less -S vcf/HanXRQChr01.vcf.gz
```

Try to find an indel. Do you see any sites with more than 1 alternate alleles? 

Pick a site and figure out, whats the minor allele frequency? How many samples were genotyped? 

We've now called variants for a single chromosome, but there are other chromosomes. In this case there are only three, but for many genomes there will be thousands of contigs. Your next challenge is to write a loop to create the genomicsdb file and then VCF for each chromosome. 

Once you have three VCF files, one for each chromosome, you can concatenate them together to make a single VCF file. We're going to use _bcftools_ which is a very fast program for manipulating vcfs as well as bcfs (the binary version of a vcf).

```bash

bcftools concat \
  vcf/HanXRQChr01.vcf.gz \
  vcf/HanXRQChr02.vcf.gz \
  vcf/HanXRQChr03.vcf.gz \
  -O z > vcf/full_genome.vcf.gz

```

You've done it! We have a VCF. Tomorrow we will fliter this file and use it for some analyses.

### Coding challenge
* Use command line tools to extract a list of all the samples in your VCF file, from the vcf file itself. They should be one name per line.
* Take the original vcf file produced and create a vcf of only high biallelic SNPs for ANN samples. 
* Use bcftools to filter your vcf file and select for sites with alternate allele frequencies > 0.01, including multi-allelic sites. 

### Daily assignments
1. Another program that is useful for filtering and formatting vcf files is [vcftools](https://vcftools.github.io/index.html). It is installed on the server. It can also do basic pop gen stats. Use it to calculate Fst between samples with ARG and ANN names.
2. You're trying to create a very stringent set of SNPs. Based on the site information GATK produces, what filters would you use? Include the actual GATK abbreviations.
3. What is strand bias and why would you filter based on it?
