# Topic 8: SNP calling with GATK
In this tutorial we're going to call SNPs with GATK. We will run steps as we talk about them 
```bash
screen -r
mkdir /home/ubuntu/log

#First we need to index the reference for GATK.
/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -jar /home/ubuntu/bin/picard.jar CreateSequenceDictionary R=/home/ubuntu/ref/HA412.v1.1.bronze.20141015.fasta O=/home/ubuntu/ref/HA412.v1.1.bronze.20141015.dict
samtools faidx HA412.v1.1.bronze.20141015.fasta

#Since we're going to apply most functions to each sample, lets make a list of samplenames
ls /home/ubuntu/fastq | grep _R1.fastq | sed s/_R1.fastq.gz//g > /home/ubuntu/samplelist.txt


#Now lets mark PCR duplicates. 
#Keep in mind, for GBS data, where the start and stop sites are set, there will be excessive false positives, so we should not mark duplicates. Here we run through it for practice but we do not use the product.
while read name \
do \
/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -jar /home/ubuntu/bin/picard.jar MarkDuplicates I=/home/ubuntu/bam/${name}.ngm.rg.clean.bam O=/home/ubuntu/bam/${name}.tmp.bam RGID=P1-1 M=/home/ubuntu/bam/${name}.duplicateinfo.txt \
done < /home/ubuntu/samplelist.txt

#We now want to create a list of potential indel sites
ls bam | grep ngm.rg.clean.bam | sed "s/^/\/home\/ubuntu\/bam\//g" > /home/ubuntu/biol525D.bamfiles.list
/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -Xmx8g -jar /home/ubuntu/bin/GenomeAnalysisTK.jar \
   -T RealignerTargetCreator \
   -R /home/ubuntu/ref/HA412.v1.1.bronze.20141015.fasta \
   -I /home/ubuntu/biol525D.bamfiles.list \
   -nt 2 \
   -log /home/ubuntu/log/biol525D.RealignerTargetCreator.log \
   -o /home/ubuntu/biol525D.realign.intervals

#Now with the list of indel sites, we need to realign in those regions
while read name \
do \
/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java -Xmx8g -jar /home/ubuntu/bin/GenomeAnalysisTK.jar \
  -T IndelRealigner \
  -R /home/ubuntu/ref/HA412.v1.1.bronze.20141015.fasta \
  -I /home/ubuntu/bam/${name}.ngm.rg.clean.bam \
  -targetIntervals /home/ubuntu/biol525D.realign.intervals \
  -o /home/ubuntu/bam/${name}.ngm.rg.clean.realign.bam \
  -log  /home/ubuntu/log/${name}.IndelRealigner.log
done < /home/ubuntu/samplelist.txt

#Next we use the haplotypecaller to create a gvcf file for each sample

#With all the gvcf files created, we jointly genotype all samples to produce a single vcf

#Lets take a look at the vcf file

#Some of those calls aren't very good. Also indels can be difficult to deal with, so lets get rid of them.

#Finally, vcf is often a difficult format to use, so lets convert it to a flat tab-separated format.
