# Topic 4: Sequence alignment

The first step is to install several programs we will be using. These should be installed in the "/home/ubuntu/bin" directory.

####Programs to install
* Samtools
* htslib
* NextGenMap
* BWA
* Picard tools

```bash
screen -r
#Download and extract samtools and htslib
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
wget https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2
tar xvjf samtools-1.3.1.tar.bz2
tar xvjf htslib-1.3.1.tar.bz2

#Install both programs
cd samtools-1.3.1 
make
sudo make install
cd ../htslib-1.3.1
make
sudo make install
cd ..

#Download and install BWA
git clone https://github.com/lh3/bwa.git
cd bwa
make
cd ..

#Download and install NextGenMap
wget https://github.com/Cibiv/NextGenMap/archive/0.5.0.tar.gz -O NextGenMap-0.5.0.tar.gz
tar xzvf NextGenMap-0.5.0.tar.gz
cd NextGenMap-0.5.0/
mkdir -p build/
cd build/
cmake ..
make
cd ../../

#Download Picardtools
wget https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar

#Index the reference for both programs
cd /home/ubuntu/ref/
/home/ubuntu/bin/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm -r HanXRQr1.0-20151230.fa
#/home/ubuntu/bin/bwa/bwa index HanXRQr1.0-20151230.fa #NOTE: This is already done because it takes an hour.
cd ..

#Test alignment on NGM
mkdir bam
/home/ubuntu/fastq/P1-1_R2.fastq.gz  -o /home/ubuntu/bam/P1-1.ngm.sam -t 2
/home/ubuntu/bin/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm -r /home/ubuntu/ref/HanXRQr1.0-20151230.fa -1 /home/ubuntu/fastq/P1-1_R1.fastq.gz -2 

###Run up to here at start of class.
#Lets examine that bam file
cd bam
less -S P1-1.ngm.sam
samtools view -bh P1-1.ngm.sam | samtools sort > P1-1.ngm.bam #Convert to bam and sort
samtools index P1-1.ngm.bam
samtools tview P1-1.ngm.bam
#Go to HanXRQChr01:2888693
#Check alignment numbers
samtools flagstat P1-1.ngm.bam > P1-1.ngm.stats.txt

#Test alignment on BWA
/home/ubuntu/bin/bwa/bwa mem /home/ubuntu/ref/HanXRQr1.0-20151230.fa /home/ubuntu/fastq/P1-1_R1.fastq.gz /home/ubuntu/fastq/P1-1_R2.fastq.gz -t 2 | samtools view -bh | samtools sort > /home/ubuntu/bam/P1-1.bwa.bam 

samtools index P1-1.bwa.bam
samtools flagstat P1-1.bwa.bam > P1-1.bwa.stats.txt

#Compare the two alignment programs.
#We'll choose NGM for this exercise.
#We then have to process the bam file for future use using picardtools

/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -jar /home/ubuntu/bin/picard.jar AddOrReplaceReadGroups I=/home/ubuntu/bam/P1-1.ngm.bam O=/home/ubuntu/bam/P1-1.ngm.rg.bam RGID=P1-1 RGLB=biol525D RGPL=ILLUMINA RGPU=biol525D RGSM=P1-1 SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0
/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -jar /home/ubuntu/bin/picard.jar CleanSam I=/home/ubuntu/bam/P1-1.ngm.rg.bam O=/home/ubuntu/bam/P1-1.ngm.rg.clean.bam VALIDATION_STRINGENCY=LENIENT

```
We now need to write a script to put these parts together and run them in one go. 
HINTS:
* Use variables for directory paths "bwa=/home/ubuntu/bin/bwa/bwa"
* Use a loop.

<details> 
  <summary>**Answer 1**  </summary>
   ```bash
   #First set up variable names
   bam=/home/ubuntu/bam
   java=/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java
   ngm=/home/ubuntu/bin/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm
   bin=/home/ubuntu/bin
   ref=/home/ubuntu/ref/HanXRQr1.0-20151230.fa
   #Then get a list of sample names, without suffixes
   ls $bam | grep ngm.rg.clean.bam | sed s/.ngm.rg.clean.bam//g > $bam/samplelist.txt
   while read name
   do
   $ngm -r $ref -1 ${name}_R1.fastq.gz -2 ${name}_R2.fastq.gz
```
</details>

	java $javarules -jar $picardtools/picard.jar SortSam I=$bam/$name.merged.bam O=$bam/$name.merged.v1.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0
	java $javarules -jar $picardtools/picard.jar AddOrReplaceReadGroups I=$bam/$name.merged.v1.bam O=$bam/$name.merged.v2.bam \
	RGID=$name RGLB=$project RGPL=ILLUMINA RGPU=$project RGSM=$name SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT COMPRESSION_LEVEL=0
	java $javarules -jar $picardtools/picard.jar CleanSam I=$bam/$name.merged.v2.bam O=$bam/$name.bam VALIDATION_STRINGENCY=LENIENT

#First we need to get a list of all the sample names


