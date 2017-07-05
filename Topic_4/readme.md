# Topic 4: Sequence alignment

The first step is to install several programs we will be using. These should be installed in the "/home/ubuntu/bin" directory.

#### Programs to install
* Samtools
* htslib
* NextGenMap
* BWA
* Picard tools

```bash
byobu

mkdir bin
cd bin

#Download and extract java, samtools and htslib
wget -c --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jre-8u131-linux-x64.tar.gz
wget https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.5/htslib-1.5.tar.bz2
tar xvjf samtools-1.5.tar.bz2
tar xvjf htslib-1.5.tar.bz2
tar xzf jre-8u131-linux-x64.tar.gz

#Install both programs
cd samtools-1.5
make
sudo make install
cd ../htslib-1.5
make
sudo make install
cd ..

#Download and install BWA
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

#Make sure cmake is installed, so that you can install NextGenMap
sudo apt-get install cmake
sudo apt-get install build-essential

#Download and install NextGenMap
wget https://github.com/Cibiv/NextGenMap/archive/v0.5.2.tar.gz -O NextGenMap-0.5.2.tar.gz
tar xzvf NextGenMap-0.5.2.tar.gz
cd NextGenMap-0.5.2/
mkdir -p build/
cd build/
cmake ..
make
cd ~

#Download Picardtools
wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar

#Index the reference for both programs. This is already done, but this is how you would do it.
#cd /home/ubuntu/ref/
#/home/ubuntu/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa
#/home/ubuntu/bin/bwa/bwa index HA412.v1.1.bronze.20141015.fasta #NOTE: This is already done because it takes an hour.
#cd ..

#Test alignment on NGM
mkdir bam
/home/ubuntu/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm \
-r /mnt/data/ref/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa \
-1 /mnt/data/Topic4/Sample1_R1.fastq.gz \
-2 /mnt/data/Topic4/Sample1_R2.fastq.gz \
-o /home/ubuntu/bam/Sample1.ngm.sam -t 2 \
--rg-id Sample1 \
--rg-sm Sample1 \ 
--rg-pl illumina \
--rg-pu Biol525D \
--rg-lb Sample1_lib

###Run up to here at start of class.
#Lets examine that bam file
cd bam
less -S Sample1.ngm.sam
samtools view -bh Sample1.ngm.sam | samtools sort > Sample1.ngm.bam #Convert to bam and sort
samtools index Sample1.ngm.bam
samtools tview Sample1.ngm.bam
#Go to groupIV:27900
#Check alignment numbers
samtools flagstat Sample1.ngm.bam > Sample1.ngm.stats.txt

#Test alignment on BWA
/home/ubuntu/bin/bwa/bwa mem /mnt/data/ref/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa \
/mnt/data/Topic4/Sample1_R1.fastq.gz \
/mnt/data/Topic4/Sample1_R2.fastq.gz \
-t 1 \
-R '@RG\tID:Sample1\tSM:Sample1\tPL:illumina\tPU:Biol525D\tLB:Sample1_lib' |\
samtools view -bh |\
samtools sort > /home/ubuntu/bam/Sample1.bwa.bam 

samtools index Sample1.bwa.bam
samtools flagstat Sample1.bwa.bam > Sample1.bwa.stats.txt

#Compare the two alignment programs.
#We'll choose NGM for this exercise.


```
We now need to write a script to put these parts together and run them in one go. 
HINTS:
* Use variables for directory paths "bwa=/home/ubuntu/bin/bwa/bwa"
* Use a loop.

<details> 
<summary> <b>Answer</b>  </summary>
  
   ```bash
   #First set up variable names
   bam=/home/ubuntu/bam
   fastq=/mnt/data/
   java=/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java
   ngm=/home/ubuntu/bin/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm
   bin=/home/ubuntu/bin
   ref=/home/ubuntu/ref/HA412.v1.1.bronze.20141015.fasta
   #Then get a list of sample names, without suffixes
   ls $fastq | grep R1.fastq | sed s/_R1.fastq.gz//g > $bam/samplelist.txt
   #Then loop through the samples
   while read name
   do
        $ngm -r $ref -1 $fastq/${name}_R1.fastq.gz -2 $fastq/${name}_R2.fastq.gz -o $bam/${name}.ngm.sam -t 2
        samtools view -bh $bam/${name}.ngm.sam | samtools sort > $bam/${name}.ngm.bam

   done < $bam/samplelist.txt
```
</details>

After your final bam files are created, and you've checked that they look good, you should remove intermediate files to save space. 
### By topic 7, you should have created cleaned bam files for all samples.

## Daily assignments
1. Is an alignment with a higher percent of mapped reads always better than one with a lower percent? Why or why not?
2. I want to reduce the percent of incorrectly mapped reads when using BWA. What setting or settings should I change in BWA?
3. What are two ways that could be used to evaluate which aligner is best?

