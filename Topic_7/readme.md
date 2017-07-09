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
project=biol525D
home=/home/ubuntu
gvcf=/home/ubuntu/gvcf



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
### Coding challenge 1:
* Take the original vcf file produced and create a vcf of only high quality indels for population 1 samples. Make sure that each indel is actually variable within population 1 samples.

### Daily assignments
1. Another program that is useful for filtering and formatting vcf files is [vcftools](https://vcftools.github.io/index.html). Successfully install this program and record the steps necessary. 
2. You're trying to create a very stringent set of SNPs. Based on the site information GATK produces, what filters would you use? Include the actual GATK abbreviations.
3. What is strand bias and why would you filter based on it?
