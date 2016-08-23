# Topic 7: SNP calling with GATK
In this tutorial we're going to call SNPs with GATK. We will run steps as we talk about them 
```bash
screen -r
mkdir /home/ubuntu/log
mkdir /home/ubuntu/gvcf
#Set up variables
ref=/home/ubuntu/ref/HA412.v1.1.bronze.20141015.fasta
java=/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java
picard=/home/ubuntu/bin/picard.jar
bam=/home/ubuntu/bam
gatk=/home/ubuntu/bin/GenomeAnalysisTK.jar
log=/home/ubuntu/log
project=biol525D
home=/home/ubuntu
fastq=/home/ubuntu/fastq
gvcf=/home/ubuntu/gvcf

#First we need to index the reference for GATK.
$java -jar $picard CreateSequenceDictionary R=$ref O=/home/ubuntu/ref/HA412.v1.1.bronze.20141015.dict
samtools faidx $ref

#Since we're going to apply most functions to each sample, lets make a list of samplenames
ls $fastq | grep _R1.fastq | sed s/_R1.fastq.gz//g > $home/samplelist.txt


#Now lets mark PCR duplicates. 
#Keep in mind, for GBS data there will be excessive false positives, so we should not mark duplicates. Here we run through it for practice but we do not use the product.
while read name 
do 
$java -jar $picard MarkDuplicates \
I=$bam/${name}.ngm.rg.clean.bam O=$bam/${name}.tmp.bam \
M=$bam/${name}.duplicateinfo.txt 
done </home/ubuntu/samplelist.txt

#We now want to create a list of potential indel sites
ls $bam | grep ngm.rg.clean.bam | grep -v bai | sed "s/^/\/home\/ubuntu\/bam\//g" > /home/ubuntu/${project}.bamfiles.list
$java -Xmx8g -jar $gatk \
   -T RealignerTargetCreator \
   -R $ref \
   -I /home/ubuntu/${project}.bamfiles.list \
   -nt 2 \
   -log $log/${project}.RealignerTargetCreator.log \
   -o /home/ubuntu/${project}.realign.intervals

#Now with the list of indel sites, we need to realign in those regions
while read name 
do 
$java -Xmx8g -jar $gatk \
  -T IndelRealigner \
  -R $ref \
  -I $bam/${name}.ngm.rg.clean.bam \
  -targetIntervals $home/${project}.realign.intervals \
  -o $bam/${name}.ngm.rg.clean.realign.bam \
  -log  $log/${name}.IndelRealigner.log
done <$home/samplelist.txt

#Next we use the haplotypecaller to create a gvcf file for each sample
time while read name 
do 
$java -Xmx8g -jar $gatk \
   -l INFO \
   -R $ref \
   -log $log/${project}.HaplotypeCaller.log \
   -T HaplotypeCaller \
   -I $bam/${name}.ngm.rg.clean.realign.bam \
   --emitRefConfidence GVCF \
   --max_alternate_alleles 1 \
   -variant_index_type LINEAR \
   -variant_index_parameter 128000 \
   -o $gvcf/${name}.GATK.g.vcf
done <$home/samplelist.txt

#Check your gvcf files to make sure each has a .idx file. If the haplotypecaller crashes, it will produce a truncated gvcf file that will eventually crash the genotypegvcf step. Note that if you give genotypegvcf a truncated file without a idx file, it will produce an idx file itself, but it still won't work. 
#The haplotypecaller tends to crash when run on multiple cores when it runs out of ram in an unpredictable fashion. 

#With all the gvcf files created, we jointly genotype all samples to produce a single vcf
ls $gvcf | grep "vcf" | grep -v ".idx"   > $home/GVCFs.samplelist.txt
tmp=""
while read prefix
do
        tmp="$tmp --variant $gvcf/$prefix"
done < $home/GVCFs.samplelist.txt

mkdir $home/tmp
$java -Xmx8g -Djava.io.tmpdir=$home/tmp -jar $gatk \
        -nt 1 \
        -l INFO \
        -R $ref \
        -log $log/GenotypeGVCFs.log \
        -T GenotypeGVCFs \
        $tmp \
        -o $home/$project.vcf \
        -hets 0.01 \
        --max_alternate_alleles 4

#Lets take a look at the vcf file
less -S $home/$project.vcf
#Lets filter this a bit. Only keeping biallelic snps and sites where less than 20% of samples are not genotyped.
#This step also lets you filter based on other quality metrics.
$java -jar $gatk \
-R $ref \
-T SelectVariants \
-V $home/$project.vcf \
-o $home/$project.snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
--maxNOCALLfraction 0.2 \
-log $log/${project}.selectvariants.log

#Finally, vcf is often a difficult format to use, so lets convert it to a flat tab-separated format.
$java -jar $gatk \
-R $ref \
-T VariantsToTable \
-V $home/$project.snps.vcf \
-F CHROM \
-F POS \
-GF GT \
-log $log/variantstotable.log \
-o $home/$project.snps.tab

#Now to filter and reformat. We want to remove the GT from the sample name, and also remove lines with *, which indicate deletions.
cat $home/$project.snps.tab |sed 's/.GT   /  /g' | sed 's/.GT$//g' | sed 's|/||g' | sed 's/\.\./NN/g' | grep -v '*' > $home/$project.snps.formatted.tab
#Note the sed 's/.GT   /  /g' command requires that the spaces are actually tabs. When copying and pasting, they are often substituted for spaces. To put an actual tab in the command, press ctrl-v, tab. 

```
###Coding challenge 1:
* Take the original vcf file produced and create a vcf of only high quality indels for population 1 samples. Make sure that each indel is actually variable within population 1 samples.

###Daily assignments
1. Another program that is useful for filtering and formatting vcf files is [vcftools](https://vcftools.github.io/index.html). Successfully install this program and record the steps necessary. 
2. You're trying to create a very stringent set of SNPs. Based on the site information GATK produces, what filters would you use? Include the actual GATK abbreviations.
3. What is strand bias and why would you filter based on it?
