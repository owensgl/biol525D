# Topic 10: Phylogenetics

The first step is to convert the vcf file into a fasta file. We are going to use PGDSpider so we need to make a spid file again, described [here](https://youtu.be/Pm6wX6IFs2I)

```
java=/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java
$java -jar bin/PGDSpider_2.1.0.3/PGDSpider2-cli.jar -inputfile biol525D.snps.vcf -inputformat VCF -outputfile biol525D.snps.fasta -outputformat FASTA -spid vcf_to_fasta.spid

#Next we should install IQ-tree
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.4.3/iqtree-omp-1.4.3-Linux.tar.gz
