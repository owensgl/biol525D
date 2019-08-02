---
title: "Topic 9 (continued): Principal Component Analysis in R"
topickey: 9.1
topictitle: "PCA"
---

Another way of visualizing the relationships between your samples is to use Principal Component Analysis. We're going to use SNPRelate for this.

Lets start fresh here so we don't have any of our previous variables interferring with this new code so open a new Rscript file and clear all objects from the workspace (the broom under the "Environment" tab on the right).

```r
#First we install SNPRelate
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("SNPRelate"))

#Then load them
library(SNPRelate)
library(tidyverse)
```
SNPRelate works with a compressed version of a genotype file called a "gds". We have to convert our vcf into a gds as the first step.

```r

snpgdsVCF2GDS("vcf/full_genome.filtered.vcf.gz",
              "vcf/full_genome.filtered.gds",
              method="biallelic.only")

#####
#Start file conversion from VCF to SNP GDS ...
#Method: exacting biallelic SNPs
#Number of samples: 10
#Parsing "vcf/full_genome.filtered.vcf.gz" ...
#	import 6750 variants.
#+ genotype   { Bit2 10x6750, 16.5K } *
#Optimize the access efficiency ...
#Clean up the fragments of GDS file:
#    open the file 'vcf/full_genome.filtered.gds' (53.1K)
#    # of fragments: 46
#    save to 'vcf/full_genome.filtered.gds.tmp'
#    rename 'vcf/full_genome.filtered.gds.tmp' (52.8K, reduced: 312B)
#    # of fragments: 20
```



```r
#Load the gds file
genofile <- snpgdsOpen("vcf/full_genome.filtered.gds")
#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)
#SNP pruning based on LD:
#Excluding 313 SNPs (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
#Working space: 10 samples, 6,437 SNPs
#    using 1 (CPU) core
#    sliding window: 500,000 basepairs, Inf SNPs
#    |LD| threshold: 0.2
#    method: composite
#Chromosome HanXRQChr01: 0.21%, 7/3,398
#Chromosome HanXRQChr02: 1.18%, 10/846
#Chromosome HanXRQChr03: 0.16%, 4/2,506
#21 markers are selected in total.
```

Holy moly, linkage is very high! In this case, its pruned ~6000 SNPs down to ~20 because they were all so closely linked. Most normal datasets will remove much fewer sites and keep much more data. In my experience, Whole Genome Sequence data can be reduced to ~20% of the original size with this procedure. 


```r
#Make a list of sites we're keeping.
snpset.id <- unlist(snpset_pruned)
#Run the PCA
pca <- snpgdsPCA(genofile, num.thread = 2, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.10, maf = 0.05,autosome.only = F)
#Principal Component Analysis (PCA) on genotypes:
#Excluding 1 SNP (monomorphic: TRUE, MAF: 0.05, missing rate: 0.1)
#Working space: 10 samples, 20 SNPs
#    using 2 (CPU) cores
#PCA:    the sum of all selected genotypes (0,1,2) = 269
#CPU capabilities: Double-Precision SSE2
#Wed Jul 24 15:54:33 2019    (internal increment: 48740)
#[==================================================] 100%, completed, 0s  
#Wed Jul 24 15:54:33 2019    Begin (eigenvalues and eigenvectors)
#Wed Jul 24 15:54:33 2019    Done.
```
PCA is done and stored in the pca object. Try printing the pca object to see whats in it.
```r
pca
```

It's got a bunch of different parts. This is useful in case we need to refer back to different elements of the results but makes it harder to plot. Lets extract some relevant parts of it.


```r
#Here's the percent variance explained for each eigenvector
pc.percent <- pca$varprop*100
round(pc.percent, 2)
# [1] 28.07 19.34 15.85  9.87  8.24  8.09  5.96  3.43  1.15  0.00
```

Now for actually plotting the eigenvectors.

```r


#Make a dataframe of your PCA results
tab <- data.frame(sample = pca$sample.id,
                  PC1 = pca$eigenvect[,1],    # the first eigenvector
                  PC2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)


#Plot a PCA image
tab %>%
  ggplot(.,aes(x=PC1,y=PC2)) + geom_point()
```

![](pca_1.jpg)

That looks cool, but isn't particularly informative because we don't actually know what each sample is. In this case the main division between the samples is ANN vs ARG, which is part of the sample name. We can leaverage this to divide our samples into groups using the _mutate()_ command. We're using _substr()_ to extract the first three characters from the sample name which will be our group ID.

```r
tab %>%
  mutate(group = substr(sample,1,3)) %>%
  ggplot(.,aes(x=PC1,y=PC2)) + 
  geom_point(aes(color=group))
```

![](pca_2.jpg)

Now we can actually tell that there's a difference between the groups. The second PC separates different samples within the ARG group, which could be useful to know. 



Plotting challenge 1
--------------------

-   Plot the 3rd and 4th principal components.



Plotting challenge 2
--------------------

-   Add individual sample labels to each point on the plot. 

HINT:
  * Try the ggrepel package
  {: .spoiler}

Lastly, lets move on to [plotting Fst across the genome](./fst.md)
