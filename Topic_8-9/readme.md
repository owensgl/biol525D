---
title: "Topics 8: Population genomics"
permalink: /Topic_8-9/
topickey: 8
topictitle: "Population genomics"
---

## Accompanying material
* [Slides](./Topic 8.pdf)

## Daily Questions:
1. For a site that is invariant in both populations (i.e. a locus with no variation), what is Fst?
2. What effect does missing data have on a PCA?
3. What is the average Fst between ARG and ANN samples in our dataset? Hint, SNPrelate can calculate Fst.

### NOTE:
* If you didn't complete creating _full_genome.vcf.gz_ in Topic 7, you can copy it to ~/vcf from /mnt/data/vcf


Last topic we called variants across the three chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~7X depth), but can be due to other factors like divergence been sample and reference. We also have indels, and SNPs with more than two alleles. Many programs strictly require biallelic sites so lets first filter the VCF to a smaller set of usable sites.
We're going to use _bcftools_ a very fast program for processing and filtering VCF files. Here are what we want to filter for:
* At least 9/10 samples genotyped (which is 18 alleles since they're diploid).
* Only one alternate allele.
* No indels.
* At least 2 copies of the alternate allele

```bash
bcftools  view \
	-c 2 \
	-i 'INFO/AN >= 18' \
	-m2 \
	-M2 \
	-v snps \
	vcf/full_genome.vcf.gz \
	-O z > vcf/full_genome.filtered.vcf.gz

#We also should index the vcf file for future use
tabix -p vcf vcf/full_genome.filtered.vcf.gz
```

##Coding challenge
* How many sites remain in the filtered VCF? How many in the chromosome HanXRQChr01?

A common first pass analysis is to use structure to look at clustering in your data. Admixture is similar to STRUCTURE but orders of magnitude faster. We're going use that, but before that we have to convert our VCF to the bed format. We're going to use plink to do that. Plink is a large set of tools for manipulating genetic data, running GWAS and calculating various stats. It's geared towards human data, so sometimes you have to force it to work with non-human data. For example, it assumes you have human chromosomes (eg 23 of them) and will complain if that doesn't see them.


```bash
cd ~/
mkdir analysis
/mnt/bin/plink --make-bed \
	--vcf vcf/full_genome.filtered.vcf.gz \
	--out vcf/full_genome.filtered \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr
```
This produces files with the suffix .nosex, .log, .fam, .bim, .bed. We can use these in Admixture.

NOTE: When using admixture you should filter your VCF for linkage (i.e. remove highly linked sites). We're going to do this later during the PCA step, so for now we're using our whole set. If you can't filter for linkage, subsetting the site also helps (i.e. selecting every 10th site).

```bash 
/mnt/bin/admixture_linux-1.3.0/admixture full_genome.filtered.bed 2
```
Uh oh that doesn't work, it produces this error message.
```bash
****                   ADMIXTURE Version 1.3.0                  ****
****                    Copyright 2008-2015                     ****
****           David Alexander, Suyash Shringarpure,            ****
****                John  Novembre, Ken Lange                   ****
****                                                            ****
****                 Please cite our paper!                     ****
****   Information at www.genetics.ucla.edu/software/admixture  ****

Random seed: 43
Point estimation method: Block relaxation algorithm
Convergence acceleration algorithm: QuasiNewton, 3 secant conditions
Point estimation will terminate when objective function delta < 0.0001
Estimation of standard errors disabled; will compute point estimates only.
Invalid chromosome code!  Use integers.
```
Our chromosomes are named HanXRQChr01, HanXRQChr02, HanXRQChr03, not integers like admixture is expecting. Like many programs, this is coded for human data where chromosomes are known and numbered clearly. We need to rename the chromosome column of the vcf so that they're integers. In this case, that means removing "HanXRQChr" from any line that starts with that, although it would depend on how your chromosomes are named.

```bash
zcat vcf/full_genome.filtered.vcf.gz |\
	sed s/^HanXRQChr//g |\
	gzip > vcf/full_genome.filtered.numericChr.vcf.gz
	
/mnt/bin/plink --make-bed \
	--vcf vcf/full_genome.filtered.numericChr.vcf.gz \
	--out vcf/full_genome.filtered.numericChr \
	--set-missing-var-ids @:# \
	--double-id \
	--allow-extra-chr

/mnt/bin/admixture_linux-1.3.0/admixture --cv vcf/full_genome.filtered.numericChr.bed 2
```
This works! With only 10 samples and ~6500 SNPs it finished almost completely. We only ran it for one value of K (2) but we should also test different K values and select the best K value.
```bash 
for K in 1 2 3 4 5; \
do /mnt/bin/admixture_linux-1.3.0/admixture --cv vcf/full_genome.filtered.numericChr.bed $K |\
tee full_genome.filtered.numericChr.${K}.out; \
done
#NOTE: "tee" takes the output of a command and saves it to a file, while 
# also printing letting it print to the screen. So we can watch the progress while also 
# saving that output. 

#Now move all the output files to the analysis directory
mv full_genome.filtered.numericChr* analysis/
```
The best K value for Admixture is typically the K value with the lowest cross-validation (CV) error. The CV error are in the .out files we saved. One easy way to look at all those scores is to print all the .out files and then keep only the lines that include "CV" using grep. 

```
cat analysis/*out | grep CV
CV error (K=1): 1.57068
CV error (K=2): 1.18264
CV error (K=3): 1.26501
CV error (K=4): 1.48708
CV error (K=5): 1.61393
```
This shows that the lowest CV error is with K=2, although K=3 is our second choice. To see how this lines up lets look at the .Q file, which shows group assignment for each sample. The Q file doesn't include sample names so we can put those together using "paste"

```bash
paste samplelist.txt analysis/full_genome.filtered.numericChr.2.Q
ANN1133	0.000010 0.999990
ANN1134	0.000010 0.999990
ANN1337	0.000010 0.999990
ANN1338	0.000010 0.999990
ANN1373	0.000010 0.999990
ARG0010	0.999990 0.000010
ARG0015	0.999990 0.000010
ARG0016	0.999990 0.000010
ARG0018	0.999990 0.000010
ARG0028	0.999990 0.000010
```
All the ANN samples are one group and all the ARG are a different group, which makes sense, since they are different species. We're going to plot these results, but before we leave the command line, lets also calculate Fst between the groups (or species in this case) using the perl tool vcf2fst.pl. This is a custom script from Greg, since we want to keep the numerator and denominator from the Fst calculations, which is hard to do.

We need two files, a sample info file and a group file. The sample info file tells the program which population each sample is in and the group file tells the program which populations to compare. We can make them here:

```
for i in `cat samplelist.txt`; 
	do echo -e "$i\t${i/%????/}"; 
done > sampleinfo.txt
echo -e "ANN\t1\nARG\t2" > popinfo.txt

zcat vcf/full_genome.filtered.vcf.gz |\
perl /mnt/bin/vcf2fst.pl sampleinfo.txt popinfo.txt \
> analysis/full_genome.filtered.fst.txt

```
We're going to move from the command line to desktop Rstudio, but as a last step lets copy our samplelist file to the analysis directory so we can use it later.
```bash
cp samplelist.txt analysis/
```

Forward to plotting on the [next page](./plotting_structure.md).

