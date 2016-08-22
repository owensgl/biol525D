# Topic 10: Phylogenetics

The first step is to convert the vcf file into a fasta file. We are going to use a perl script in this github directory. Download SNPtable2fasta.pl and move it into your bin directory. It converts the tab style file that we created to a fasta, for heterozygous sites it picks one random allele. IQtree and Splitstree both have seem to ignore ambiguous bases, so we're not going to use them.

We also have to remove invariant sites, so that we can use an ascertainment bias in IQtree. 
IQtree also has a [webserver](http://iqtree.cibiv.univie.ac.at/).

```
cat biol525D.snps.formatted.tab | perl bin/SNPtable2fasta_pickoneallele.pl | perl bin/fasta2removeinvariant.pl > biol525D.snps.fasta

#Next we should install IQ-tree
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.4.3/iqtree-omp-1.4.3-Linux.tar.gz

#Next we run it
bin/iqtree-omp-1.4.3-Linux/bin/iqtree-omp -s biol525D.snps.fasta -st DNA -m TEST+ASC -nt 2 -alrt 1000 -lbp 1000

#This produces several output files, including a log and a couple different versions of the treefile. The -m TEST command does a model test and selects the best substition model. 


