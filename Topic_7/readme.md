# Topic 7: SNP calling with GATK
In this tutorial we're going to call SNPs with GATK. We will run steps as we talk about them 
```bash
byobu
mkdir /home/ubuntu/log
mkdir /home/ubuntu/gvcf
#Set up variables
ref=/mnt/data/ref/Gasterosteus_aculeatus.BROADS1.dna_rm.toplevel.fa
java=/home/ubuntu/bin/jre1.8.0_131/bin/java
bam=/home/ubuntu/bam
gatk=/mnt/data/bin/GenomeAnalysisTK.jar
log=/home/ubuntu/log
home=/home/ubuntu
gvcf=/home/ubuntu/gvcf
fastq=/mnt/data/Topic4
ngm=/home/ubuntu/bin/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm
bin=/home/ubuntu/bin
project=Biol525D
#We need to make aligned bam files for all our samples. I've reduced the data in each sample because of RAM limitations, so please rerun this command
for name in `ls $fastq | grep R1.fastq.gz | sed s/_R1.fastq.gz//g`
do
     $ngm \
       -r $ref \
       -1 $fastq/${name}_R1.fastq.gz \
       -2 $fastq/${name}_R2.fastq.gz \
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
     rm $bam/${name}.ngm.sam

done

#Since we're going to apply most functions to each sample, lets make a list of samplenames
ls $fastq | grep _R1.fastq | sed s/_R1.fastq.gz//g > $home/samplelist.txt


#Now lets mark PCR duplicates. 
#Keep in mind, for GBS data there will be excessive false positives, so we should not mark duplicates. Here we run through it for practice but we do not use the product.
while read name 
do 
$java -jar $picard MarkDuplicates \
I=$bam/${name}.ngm.rg.clean.bam O=$bam/${name}.tmp.bam \
M=$log/${name}.duplicateinfo.txt 
done < $home/samplelist.txt



#Next we use the haplotypecaller to create a gvcf file for each sample
while read name 
do 
$java -Xmx3g -jar $gatk \
   -l INFO \
   -R $ref \
   -log $log/${project}.HaplotypeCaller.log \
   -T HaplotypeCaller \
   -I $bam/${name}.ngm.bam \
   --emitRefConfidence GVCF \
   --max_alternate_alleles 2 \
   -o $gvcf/${name}.GATK.g.vcf
done <$home/samplelist.txt

#Check your gvcf files to make sure each has a .idx file. If the haplotypecaller crashes, it will produce a truncated gvcf file that will eventually crash the genotypegvcf step. Note that if you give genotypegvcf a truncated file without a idx file, it will produce an idx file itself, but it still won't work. 
#The haplotypecaller tends to crash when run on multiple cores when it runs out of ram in an unpredictable fashion. 

#With all the gvcf files created, we jointly genotype all samples to produce a single vcf
ls $gvcf | grep "vcf" | grep -v ".idx"   > $home/GVCFs.samplelist.txt

#This loop is going to make a string that has all the sample names, so we can give that to GATK.
tmp=""
while read prefix
do
        tmp="$tmp --variant $gvcf/$prefix"
done < $home/GVCFs.samplelist.txt

mkdir $home/tmp
$java -Xmx3g -Djava.io.tmpdir=$home -jar $gatk \
        -nt 1 \
        -l INFO \
        -R $ref \
        -log $log/GenotypeGVCFs.log \
        -T GenotypeGVCFs \
        $tmp \
        -o $home/$project.vcf \
        -hets 0.01 \
        --max_alternate_alleles 4


####RUN UP TO THIS STEP AT THE START OF THE CLASS

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
--maxNOCALLfraction 0.4 \
-log $log/$project.selectvariants.log

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

#We don't have very many reads, so for many of the genotypes we can't confidently call SNPs. One way to deal with this is to take into account the uncertainty of the genotype call using ANGSD.
#First we install ANGSD
cd ~/bin
git clone https://github.com/ANGSD/angsd.git 
cd angsd ;make HTSSRC=../htslib-1.5

#ANGSD has many different functions, but for now lets just calculate allele frequencies for our three samples
ls -d $PWD $bam/*.* | grep .bam$ | grep -v bai > ~/$project.bamlist.txt
~/bin/angsd/angsd -out $project.angsd -doMajorMinor 2 -doMaf 8 -bam  ~/$project.bamlist.txt -doCounts 1  -minMaf 0.1

#This is much faster than GATK, although it does not do local realignment of haplotypes. 



```
### Coding challenge
* Take the original vcf file produced and create a vcf of only high quality indels for population 1 samples. Make sure that each indel is actually variable within population 1 samples.
* It can be better to use genotype likelihoods instead of called genotypes when you have low depth. One tool for doing this is [ngsTools](https://github.com/mfumagalli/ngsTools#ngscovar). Install this program and run a PCA on your samples.

### Daily assignments
1. Another program that is useful for filtering and formatting vcf files is [vcftools](https://vcftools.github.io/index.html). Successfully install this program and record the steps necessary. 
2. You're trying to create a very stringent set of SNPs. Based on the site information GATK produces, what filters would you use? Include the actual GATK abbreviations.
3. What is strand bias and why would you filter based on it?
