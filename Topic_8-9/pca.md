# Principal Component Analysis in R

First used cyberduck to download "Biol525D.snps.vcf" from your server
```r
#First we install some packages
source("http://bioconductor.org/biocLite.R")
biocLite("SNPRelate")
```

```
## 
## The downloaded binary packages are in
## 	/var/folders/xr/sh3jfkts64q_lxbhbsftgmhw0000gn/T//RtmptFKnPA/downloaded_packages
```

```r
#Then load them
library(SNPRelate)
library(tidyverse)
```


```r
#Set up file names
vcf_filename<- c("Downloads/Biol525D.snps.vcf")
gds_filename<- c("Downloads/Biol525D.snps.gds")
sampleinfo_filename <- c("Downloads/Biol525D_popinfo.csv")
#Convert your vcf to gds for use with snprelate
snpgdsVCF2GDS(vcf_filename, gds_filename,  method="biallelic.only",ignore.chr.prefix="group")
```

```
## VCF Format ==> SNP GDS Format
## Method: exacting biallelic SNPs
## Number of samples: 3
## Parsing "/Users/gregoryowens/Downloads/Biol525D.snps.vcf" ...
## 	import 15705 variants.
## + genotype   { Bit2 3x15705, 11.5K } *
## Optimize the access efficiency ...
## Clean up the fragments of GDS file:
##     open the file '/Users/gregoryowens/Downloads/Biol525D.snps.gds' (128.5K)
##     # of fragments: 42
##     save to '/Users/gregoryowens/Downloads/Biol525D.snps.gds.tmp'
##     rename '/Users/gregoryowens/Downloads/Biol525D.snps.gds.tmp' (128.3K, reduced: 264B)
##     # of fragments: 20
```


```r
#Load the gds file
genofile <- snpgdsOpen(gds_filename)
#Prune for linkage
snpset_pruned <- snpgdsLDpruning(genofile, autosome.only=F)
```

```
## SNP pruning based on LD:
## Excluding 4,581 SNPs (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 3 samples, 11,124 SNPs
##     using 1 (CPU) core
## 	Sliding window: 500000 basepairs, Inf SNPs
## 	|LD| threshold: 0.2
## Chromosome IV: 4.78%, 59/1234
## Chromosome I: 5.18%, 56/1081
## Chromosome VII: 4.41%, 53/1201
## Chromosome II: 5.04%, 40/793
## Chromosome IX: 4.84%, 39/806
## Chromosome XIX: 6.39%, 34/532
## Chromosome XIII: 5.27%, 34/645
## Chromosome XX: 4.59%, 37/806
## Chromosome VIII: 4.95%, 37/748
## Chromosome XII: 5.38%, 35/651
## Chromosome XVI: 5.78%, 32/554
## Chromosome VI: 5.75%, 34/591
## Chromosome III: 5.42%, 30/554
## Chromosome XI: 4.75%, 31/652
## Chromosome XVIII: 5.05%, 30/594
## Chromosome XV: 6.85%, 33/482
## Chromosome X: 5.35%, 29/542
## Chromosome XIV: 5.13%, 28/546
## Chromosome XVII: 4.67%, 28/600
## Chromosome V: 5.14%, 22/428
## Chromosome XXI: 4.56%, 22/482
## Chromosome scaffold_27: 6.25%, 9/144
## Chromosome scaffold_37: 5.63%, 4/71
## Chromosome scaffold_47: 5.41%, 2/37
## Chromosome scaffold_48: 5.41%, 2/37
## Chromosome scaffold_54: 3.57%, 1/28
## Chromosome scaffold_56: 15.38%, 2/13
## Chromosome scaffold_58: 16.67%, 1/6
## Chromosome scaffold_61: 10.00%, 1/10
## Chromosome scaffold_67: 4.55%, 2/44
## Chromosome scaffold_68: 5.71%, 2/35
## Chromosome scaffold_69: 11.11%, 1/9
## Chromosome scaffold_74: 4.44%, 2/45
## Chromosome scaffold_76: 100.00%, 1/1
## Chromosome scaffold_80: 20.00%, 1/5
## Chromosome scaffold_84: 8.33%, 1/12
## Chromosome scaffold_89: 11.11%, 1/9
## Chromosome scaffold_90: 12.50%, 1/8
## Chromosome scaffold_88: 16.67%, 1/6
## Chromosome scaffold_95: 6.25%, 1/16
## Chromosome scaffold_99: 1.82%, 1/55
## Chromosome scaffold_98: 12.50%, 1/8
## Chromosome scaffold_121: 25.00%, 1/4
## Chromosome scaffold_101: 16.67%, 1/6
## Chromosome scaffold_106: 50.00%, 1/2
## Chromosome scaffold_114: 7.14%, 1/14
## Chromosome scaffold_112: 5.26%, 1/19
## Chromosome scaffold_115: 5.26%, 1/19
## Chromosome scaffold_111: 11.11%, 1/9
## Chromosome scaffold_122: 3.45%, 1/29
## Chromosome scaffold_120: 8.33%, 1/12
## Chromosome scaffold_126: 7.14%, 1/14
## Chromosome scaffold_129: 40.00%, 2/5
## Chromosome scaffold_133: 25.00%, 1/4
## Chromosome scaffold_152: 25.00%, 1/4
## Chromosome scaffold_132: 16.67%, 1/6
## Chromosome scaffold_128: 25.00%, 1/4
## Chromosome scaffold_137: 20.00%, 1/5
## Chromosome scaffold_135: 20.00%, 2/10
## Chromosome scaffold_130: 33.33%, 1/3
## Chromosome scaffold_146: 50.00%, 1/2
## Chromosome scaffold_139: 16.67%, 1/6
## Chromosome scaffold_148: 50.00%, 1/2
## Chromosome scaffold_151: 6.25%, 1/16
## Chromosome scaffold_150: 16.67%, 1/6
## Chromosome scaffold_149: 7.69%, 1/13
## Chromosome scaffold_175: 20.00%, 1/5
## Chromosome scaffold_157: 11.11%, 1/9
## Chromosome scaffold_156: 16.67%, 1/6
## Chromosome scaffold_165: 8.33%, 1/12
## Chromosome scaffold_159: 16.67%, 1/6
## Chromosome scaffold_163: 20.00%, 1/5
## Chromosome scaffold_161: 20.00%, 1/5
## Chromosome scaffold_168: 5.00%, 1/20
## Chromosome scaffold_169: 11.11%, 1/9
## Chromosome scaffold_177: 16.67%, 1/6
## Chromosome scaffold_184: 5.26%, 1/19
## Chromosome scaffold_195: 100.00%, 1/1
## Chromosome scaffold_201: 33.33%, 1/3
## Chromosome scaffold_180: 14.29%, 1/7
## Chromosome scaffold_181: 16.67%, 1/6
## Chromosome scaffold_178: 50.00%, 1/2
## Chromosome scaffold_182: 25.00%, 1/4
## Chromosome scaffold_197: 16.67%, 1/6
## Chromosome scaffold_202: 14.29%, 1/7
## Chromosome scaffold_208: 25.00%, 1/4
## Chromosome scaffold_206: 16.67%, 1/6
## Chromosome scaffold_199: 33.33%, 1/3
## Chromosome scaffold_216: 20.00%, 1/5
## Chromosome scaffold_204: 25.00%, 1/4
## Chromosome scaffold_225: 25.00%, 1/4
## Chromosome scaffold_218: 33.33%, 1/3
## Chromosome scaffold_214: 8.33%, 1/12
## Chromosome scaffold_213: 14.29%, 1/7
## Chromosome scaffold_211: 22.22%, 2/9
## Chromosome scaffold_219: 50.00%, 1/2
## Chromosome scaffold_243: 50.00%, 1/2
## Chromosome scaffold_281: 20.00%, 1/5
## Chromosome scaffold_240: 50.00%, 1/2
## Chromosome scaffold_229: 33.33%, 1/3
## Chromosome scaffold_236: 100.00%, 1/1
## Chromosome scaffold_238: 20.00%, 1/5
## Chromosome scaffold_249: 25.00%, 1/4
## Chromosome scaffold_246: 20.00%, 1/5
## Chromosome scaffold_256: 11.11%, 1/9
## Chromosome scaffold_300: 14.29%, 1/7
## Chromosome scaffold_326: 100.00%, 1/1
## Chromosome scaffold_323: 100.00%, 1/1
## Chromosome scaffold_338: 33.33%, 1/3
## Chromosome scaffold_274: 100.00%, 1/1
## Chromosome scaffold_270: 20.00%, 1/5
## Chromosome scaffold_273: 20.00%, 1/5
## Chromosome scaffold_277: 100.00%, 1/1
## Chromosome scaffold_349: 100.00%, 1/1
## Chromosome scaffold_354: 14.29%, 1/7
## Chromosome scaffold_382: 25.00%, 1/4
## Chromosome scaffold_363: 16.67%, 1/6
## Chromosome scaffold_364: 100.00%, 1/1
## Chromosome scaffold_381: 50.00%, 1/2
## Chromosome scaffold_425: 50.00%, 1/2
## Chromosome scaffold_440: 33.33%, 1/3
## Chromosome scaffold_478: 50.00%, 1/2
## Chromosome scaffold_475: 50.00%, 1/2
## Chromosome scaffold_522: 100.00%, 1/1
## Chromosome scaffold_524: 100.00%, 1/1
## Chromosome scaffold_542: 100.00%, 1/1
## Chromosome scaffold_560: 33.33%, 1/3
## Chromosome scaffold_604: 100.00%, 1/1
## Chromosome scaffold_621: 100.00%, 1/1
## Chromosome scaffold_707: 33.33%, 1/3
## Chromosome scaffold_719: 100.00%, 1/1
## Chromosome scaffold_718: 25.00%, 1/4
## Chromosome scaffold_869: 50.00%, 1/2
## Chromosome scaffold_740: 20.00%, 1/5
## Chromosome scaffold_757: 50.00%, 1/2
## Chromosome scaffold_762: 50.00%, 1/2
## Chromosome scaffold_891: 50.00%, 1/2
## Chromosome scaffold_926: 25.00%, 1/4
## Chromosome scaffold_1043: 100.00%, 1/1
## Chromosome scaffold_1050: 50.00%, 1/2
## Chromosome scaffold_1075: 25.00%, 1/4
## Chromosome scaffold_1250: 100.00%, 1/1
## Chromosome scaffold_1278: 33.33%, 1/3
## Chromosome scaffold_1282: 50.00%, 1/2
## Chromosome scaffold_1311: 100.00%, 1/1
## Chromosome scaffold_1331: 100.00%, 1/1
## Chromosome scaffold_1358: 33.33%, 1/3
## Chromosome scaffold_1395: 33.33%, 1/3
## Chromosome scaffold_1695: 50.00%, 1/2
## 891 SNPs are selected in total.
```

```r
snpset.id <- unlist(snpset_pruned)
#Run the PCA
pca <- snpgdsPCA(genofile, num.thread = 2, eigen.cnt = 16, snp.id = snpset.id, missing.rate = 0.10, maf = 0.05,autosome.only = F)
```

```
## Principal Component Analysis (PCA) on genotypes:
## Excluding 446 SNPs (monomorphic: TRUE, < MAF: 0.05, or > missing rate: 0.1)
## Working space: 3 samples, 445 SNPs
##     using 2 (CPU) cores
## PCA:	the sum of all selected genotypes (0, 1 and 2) = 1585
## Tue Jul 11 16:33:18 2017    (internal increment: 162472)
## 
[..................................................]  0%, ETC: ---    
[==================================================] 100%, completed      
## Tue Jul 11 16:33:18 2017    Begin (eigenvalues and eigenvectors)
## Tue Jul 11 16:33:18 2017    Done.
```

```r
#Lets take a look at the percent variance explained
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
```

```
## [1] 56.16 43.84  0.00
```


```r
#Load your sample information for plotting purposes.
sampleinfo <- read.csv(sampleinfo_filename,header=T)

#Make a dataframe of your PCA results
tab <- data.frame(name = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#Merge the sampleinfo into that
tab <- merge(tab, sampleinfo)

#Plot a PCA image
ggplot(data=tab,aes(EV1,EV2)) + geom_point()
```

![](figures/pca1-1.png)


```r
#Next lets color code by population and add axis labels
ggplot(data=tab,aes(EV1,EV2)) + geom_point(aes(color=as.factor(pop))) + ylab("Principal component 2") + xlab("Principal component 1")
```

![](figures/pca1-2.png)

```r
#We can make that look nicer
ggplot(data=tab,aes(EV1,EV2)) + geom_point(aes(color=as.factor(pop))) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + scale_color_discrete(name="Population") +
  theme(panel.border = element_rect(fill = NA, colour = "grey50")) 
```

![](figures/pca1-3.png)


Plotting challenge 1
--------------------

-   Plot the 3rd and 4th principal components and color code the points by the color sample info.



Plotting challenge 2
--------------------

-   Use the latitude and longitude information to plot your pca scores on a map. There is many ways of doing this, look it up.

Lastly, lets move on to [plotting Fst across the genome](https://github.com/owensgl/biol525D/blob/master/Topic_8-9/fst.md)
