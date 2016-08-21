# Topic 10: Phylogenetics

The first step is to convert the vcf file into a fasta file. We are going to use a perl script in this github directory. Download SNPtable2fasta.pl and move it into your bin directory. It converts the tab style file that we created to a fasta, coding heterozygous sites with IUPAC single letters.


```
cat biol525D.snps.formatted.tab | perl bin/SNPtable2fasta.pl > biol525D.snps.fa

#Next we should install IQ-tree
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.4.3/iqtree-omp-1.4.3-Linux.tar.gz
