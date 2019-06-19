---
title: Topic 7
permalink: /Topic_7/
---

# Topic 7: SNP calling with GATK
In this tutorial we're going to call SNPs with GATK. We will run steps as we talk about them 
```bash
byobu

#To make this easier, we're going to set up our directories and make a bunch of variables with paths to programs we're going to use and directories we want to make files in. 
mkdir /mnt/<USERNAME>/log
mkdir /mnt/<USERNAME>/gvcf
#Set up variables
ref=/home/biol525d/ref
java=/home/biol525d/bin/jre1.8.0_131/bin/java
bam=/mnt/<USERNAME>/bam
gatk=/home/biol525d/bin/gatk-4.0.6.0/gatk
picard=/home/biol525d/bin/picard.jar
log=/mnt/<USERNAME>/log
home=/mnt/<USERNAME>
gvcf=/mnt/<USERNAME>/gvcf
fastq=/home/biol525d/Topic_4/fastq
bin=/home/biol525d/bin
project=biol525d


#Since we're going to apply most functions to each sample, lets make a list of samplenames
ls $bam | grep .ngm.bam$ | sed s/.ngm.bam//g > $home/samplelist.txt


#Now lets mark PCR duplicates. 
#If you were using GBS data, you wouldn't want to mark duplicates.
while read name 
do 
$java -jar $picard MarkDuplicates \
I=$bam/${name}.ngm.bam O=$bam/${name}.ngm.dedup.bam \
M=$log/${name}.duplicateinfo.txt
samtools index $name.ngm.dedup.bam
done < $home/samplelist.txt

$java -jar $picard CreateSequenceDictionary R= $ref/reference.fa O= $ref/reference.dict
samtools faidx $ref/reference.fa 
#Next we use the haplotypecaller to create a gvcf file. This takes about 20 minutes per sample, so for right now we're only going to run it on one sample. 
for name in `cat $home/samplelist.txt | head -n 1`
do 
$gatk --java-options "-Xmx4g" HaplotypeCaller \
   -R $ref \
   -I $bam/${name}.ngm.dedup.bam \
   --native-pair-hmm-threads 1 \
   -ERC GVCF \
   -O $gvcf/${name}.ngm.dedup.g.vcf
done 

#Check your gvcf file to make sure it has a .idx index file. If the haplotypecaller crashes, it will produce a truncated gvcf file that will eventually crash the genotypegvcf step. Note that if you give genotypegvcf a truncated file without a idx file, it will produce an idx file itself, but it still won't work. 

#We don't have time to run each of those, so if the first file worked, copy the completed gvcf files into your directory.
cp /home/biol525d/Topic_7/gvcf/* /mnt/<USERNAME>/gvcf/

```

The next step is to import our gvcf files into a genomicsDB file. This is a compressed database representation of all the read data in our samples. It has two important features to remember:
1) Each time you call GenomicsDBImport, you create a database for a single interval. This means that you can parallelize it easier, for example by calling it once per chromosome.
2) The GenomicsDB file contains all the information of your GVCF files, but can't be added to, and can't be back transformed into a gvcf. That means if you get more samples, you can't just add them to your genomicdDB file, you have to go back to the gvcf files.


We need to create a map file to GATK where our gvcf files are and what sample is in each. Because we use a regular naming scheme for our samples, we can create that using a bash script.
This is what we're looking for:
sample1      sample1.vcf.gz
sample2      sample2.vcf.gz
sample3      sample3.vcf.gz

```bash

for i in `ls $gvcf | grep "vcf" | grep -v ".idx" | sed s/.ngm.dedup.g.vcf//g`
do
echo -e "$i\t$gvcf/$i.ngm.dedup.g.vcf"
done > $home/$project.sample_map

#Break down the piped grep and sed commands above and figure out what they do. Try removing one part and see what you.

#Next we call GenomicsDBImport to actually create the database. This will take about 10 minutes.

$gatk --java-options "-Xmx4g -Xms4g" \
       GenomicsDBImport \
       --genomicsdb-workspace-path $home/db \
       --batch-size 50 \
       -L Chr1 \
       --sample-name-map $home/$project.sample_map \
       --reader-threads 1

#With the genomicsDB created, we're finally ready to actually call variants and output a vcf. This will take around 6 minutes.

$gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $ref/reference.fa \
   -V gendb://$home/db \
   -O $home/$project.vcf.gz




#Lets take a look at the vcf file
less -S $home/$project.vcf.gz
#Try to find an indel 

#Lets filter this a bit. Only keeping biallelic snps, with a minor allele frequency > 10%. Also, sometimes you'll want to work will less SNPs than your full dataset, so lets subset to a random 10% of the total sites. 
#This step also lets you filter based on other quality metrics.
$gatk SelectVariants \
-R $ref/reference.fa \
-V $home/$project.vcf.gz \
-O $home/$project.snps.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-fraction 0.1 \
-select "AF > 0.1" 

#Finally, vcf is often a difficult format to use, so lets convert it to a flat tab-separated format.
$gatk VariantsToTable \
-R $ref/reference.fa \
-V $home/$project.snps.vcf.gz \
-F CHROM \
-F POS \
-GF GT \
-O $home/$project.snps.tab

#Now to filter and reformat. We want to remove the GT from the sample name, and also remove lines with *, which indicate deletions.
cat $home/$project.snps.tab |sed 's/.GT   /  /g' | sed 's/.GT$//g' | sed 's|/||g' | sed 's/\.\./NN/g' | grep -v '*' > $home/$project.snps.formatted.tab
#Note the sed 's/.GT   /  /g' command requires that the spaces are actually tabs. When copying and pasting, they are often substituted for spaces. To put an actual tab in the command, press ctrl-v, tab. 

#This tab delimited file is easier to parse if you want to write your own scripts.

```
### Coding challenge
* Use command line tools to extract a list of all the samples in your VCF file, from the vcf file itself. They should be one name per line.
* Take the original vcf file produced and create a vcf of only high quality indels for samples 1-25,50-75. Make sure that each indel is actually variable in those samples.
* Use GATK to filter your vcf file and select for sites with alternate allele frequencies > 0.01, including multi-allelic sites. 

### Daily assignments
1. Another program that is useful for filtering and formatting vcf files is [vcftools](https://vcftools.github.io/index.html). It is installed on the server. It can also do basic pop gen stats. Use it to calculate Fst between samples 1-50 and 51-100.
2. You're trying to create a very stringent set of SNPs. Based on the site information GATK produces, what filters would you use? Include the actual GATK abbreviations.
3. What is strand bias and why would you filter based on it?
