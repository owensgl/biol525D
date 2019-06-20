---
title: "Topic 9 (final): Calculating and plotting Fst"
topickey: 9.2
topictitle: "Plotting Fst"
---

Now to calculate Fst, we're going back to the server.

```bash


#Next we use it to calculate Fst
#We need a list of samples for the two different populations. 
cd ~
seq -f "%03g" 1 50 > samplelist.pop1.txt
seq -f "%03g" 51 100 > samplelist.pop2.txt
vcftools  \
--gzvcf biol525d.snps.vcf.gz \
--weir-fst-pop samplelist.pop1.txt --weir-fst-pop samplelist.pop2.txt \
--out biol525d
#Take a look at the fst file, does it have reasonable data? Is it all missing data?
```
**Question 1**:
* Can you calculate overall Fst from this file and if so how would you do it using command line tools?

<details> 
<summary><b>Answer 1</b>  </summary>

  
    Fst is a ratio so calculating the overall values requires summing the numerator and denominator for each locus, which we don't have. 
    
</details>

---

Next we'll take the Fst values and plot them in R. Transfer the biol525D.weir.fst file to your computer


```r
install.packages("qqman", repos = "http://cran.us.r-project.org")
library(qqman)
library(tidyverse)
```

```r
fst.filename <- "Downloads/biol525D.weir.fst"
data <- read.delim(fst.filename, header=T)
#Now one problem with plotting this is that the chromosomes are not intergers
summary(data$CHROM)
```

```
## Chr1 
## 5119 
```

```r
# This strips the "group" from the chromosome name
data$CHROM <- gsub("Chr", "", data$CHROM)

# This converts it to numeric
data$CHROM <- as.numeric(data$CHROM)

# Lets also remove values that are NA
data %>% filter(WEIR_AND_COCKERHAM_FST != "NaN") -> data

# It's important to make sure that the chromosomes are numeric instead of character
data$CHROM <- as.numeric(data$CHROM)

# Lets do a basic plot using the manhattan tool in qqman. This is generally designed
# for plotting pvalues from GWAS, but it works here.
manhattan(data, chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST", snp="POS",
          logp=FALSE,ylab="Fst")
```

![](figure/fst1-1.png)


**Question 2**:
* Earlier, we filtered low frequency variants out of the vcf file. What effect do those sites have on overall Fst?

**Coding challenges**:
* Do a manhattan plot in ggplot or base R. As a bonus, make it look nice.
* Calculate the \*mean Fst in for 5kb windows across the genome. In this case, just mean the individual Fst values, although for real data make sure you sum the numerator and denominator of Fst for windows. HINT: dplyr's group_by command can help.
* Plot the windowed Fst.




