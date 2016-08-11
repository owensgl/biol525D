# Topic 4: Sequence alignment

The first step is to install several programs we will be using. These should be installed in the "/home/ubuntu/bin" directory.

####Programs to install
* Samtools
* htslib
* NextGenMap
* BWA

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

#Index the reference for both programs
cd /home/ubuntu/ref/
/home/ubuntu/bin/NextGenMap-0.5.0/bin/ngm-0.5.0/ngm -r HanXRQr1.0-20151230.fa
/home/ubuntu/bin/bwa/bwa index HanXRQr1.0-20151230.fa
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




