---
title: "Topic 4: Sequence Alignment"
permalink: /Topic_4/
---

### Accompanying material

* [Slides](./Topic 4.pdf)


```bash
byobu

#Navigate to your working directory
cd /mnt/<USERNAME>

#Copy the reference directory to your working directory
cp -r /home/biol525d/ref/ ./

#Index the reference for BWA and NGM. This step is fast for a very small reference that we're working with but can take an hour with full genomes. 
/home/biol525d/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r ref/reference.fa
/home/biol525d/bin/bwa/bwa index ref/reference.fa


#Test alignment on NGM
mkdir bam
/home/biol525d/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm \
  -r /mnt/<USERNAME>/ref/reference.fa \
  -1 /home/biol525d/Topic_4/fastq/001.R1.fastq.gz \
  -2 /home/biol525d/Topic_4/fastq/001.R2.fastq.gz \
  -o /mnt/<USERNAME>/bam/001.ngm.sam \
  -t 2 \
  --rg-id Sample_001 \
  --rg-sm Sample_001 \
  --rg-pl illumina \
  --rg-pu biol525d \
  --rg-lb Sample_001_lib

#Lets examine that bam file
cd bam
less -S 001.ngm.sam
#Notice the @PG line that includes the program call that created the sam file. This is useful for record keeping.

```
Lets examine the sam file. It contains all the information on the reads from the fastq file, but also alignment information. 
### Questions:
1. How are reads ordered in the sam file? 
2. What does the 6th column represent? What would the string "1S93M6S" mean?
3. What are three possible reasons why mapping quality could be low for a particular read?

```bash

#The next step is to convert our same file (human readable) to a bam file (machine readable) and sort reads by their aligned position.
samtools view -bh 001.ngm.sam | samtools sort > 001.ngm.bam 
#Look at the file sizes, which is smaller? 

#Now we want to take a look at our aligned reads. First we index the file, then we use samtools tview.
samtools index 001.ngm.bam 
samtools tview 001.ngm.bam --reference /mnt/<USERNAME>/ref/reference.fa
#use ? to open the help menu. Scroll left and right with H and L. 
#Try to find positions where the sample doesn't have the reference allele. 


#Next we want to try a different aligner. This time we're going to directly pipe our output between programs instead of writing intermediate files
/home/biol525d/bin/bwa/bwa mem \
  /mnt/<USERNAME>/ref/reference.fa \
  /home/biol525d/Topic_4/fastq/001.R1.fastq.gz \
  /home/biol525d/Topic_4/fastq/001.R2.fastq.gz \
  -t 2 \
  -R '@RG\tID:Sample_001\tSM:Sample_001\tPL:illumina\tPU:biol525d\tLB:Sample_001_lib' |\
 samtools view -bh |\
 samtools sort > 001.bwa.bam 

samtools index 001.bwa.bam

#We can use flagstat to look at general stats about the alignment including how many reads aligned. This can help you pick an alignment program.
samtools flagstat 001.ngm.bam > 001.ngm.stats.txt
samtools flagstat 001.bwa.bam > 001.bwa.stats.txt


```
We have 100 samples, so we don't want to have to type these commands out 100 times. Write a bash script to produced a sorted bam file for each sample.
HINTS:
* Use variables for directory paths "bwa=/home/biol525d/bin/bwa/bwa"
* Use a loop.

<details><summary><b>Answer</b></summary><p>

```bash
    #First set up variable names
    bam=/mnt/<USERNAME>/bam
    fastq=/home/biol525d/Topic_4/fastq
   ngm=/home/biol525d/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm
   ref=/mnt/<USERNAME>/ref/reference.fa
   project=biol525d
    #Then get a list of sample names, without suffixes
    ls $fastq | grep R1.fastq.gz | sed s/.R1.fastq.gz//g > $bam/samplelist.txt
    #Then loop through the samples
    while read name
    do
         $ngm \
           -r $ref \
           -1 $fastq/${name}.R1.fastq.gz \
           -2 $fastq/${name}.R2.fastq.gz \
           -o $bam/${name}.ngm.sam \
           --rg-id $name \
           --rg-sm $name \
           --rg-pl illumina \
           --rg-pu $project \
           --rg-lb ${name}_lib \
           -t 1 
         samtools view -bh $bam/${name}.ngm.sam |\
         samtools sort > $bam/${name}.ngm.bam
         samtools index $bam/${name}.ngm.bam

    done < $bam/samplelist.txt
```
</p></details>

After your final bam files are created, and you've checked that they look good, you should remove intermediate files to save space. You can build file removal into your bash scripts, but it is often helpful to only add that in once the script works. It's hard to troubleshoot a failed script if it deletes everything as it goes. 
### By topic 7, you should have created cleaned bam files for all samples.

## Daily assignments
1. Is an alignment with a higher percent of mapped reads always better than one with a lower percent? Why or why not?
2. I want to reduce the percent of incorrectly mapped reads when using BWA. What setting or settings should I change in BWA?
3. What are two ways that could be used to evaluate which aligner is best?

