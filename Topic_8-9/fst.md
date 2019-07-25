---
title: "Topic 9 (final): Calculating and plotting Fst"
topickey: 9.2
topictitle: "Plotting Fst"
---

We're now moving onto plotting Fst values, so you should again start a new Rscript and clear your environment.


We have a the Fst values for each site in the genome. Lets load that into R.
```r
library(tidyverse)

fst <- read_tsv("analysis/full_genome.filtered.weir.fst")
fst
# A tibble: 6,750 x 3
   CHROM         POS WEIR_AND_COCKERHAM_FST
   <chr>       <dbl>                  <dbl>
 1 HanXRQChr01  1134                 1     
 2 HanXRQChr01  1137                 0.214 
 3 HanXRQChr01  1146               NaN     
 4 HanXRQChr01  1160                 0.0543
 5 HanXRQChr01  1166                 1     
 6 HanXRQChr01  1335                 0.352 
 7 HanXRQChr01  1404                 1     
 8 HanXRQChr01  1429                 0.352 
 9 HanXRQChr01  1450                 1     
10 HanXRQChr01  1458                 1     
# … with 6,740 more rows
```
That worked, but its kind of annoying to have to type in such long ALL CAPS names, so we should rename them using the _rename()_ command.
```r
fst <- fst %>%
  rename(chr=CHROM,pos=POS, fst=WEIR_AND_COCKERHAM_FST) 
fst
# A tibble: 6,750 x 3
   chr           pos      fst
   <chr>       <dbl>    <dbl>
 1 HanXRQChr01  1134   1     
 2 HanXRQChr01  1137   0.214 
 3 HanXRQChr01  1146 NaN     
 4 HanXRQChr01  1160   0.0543
 5 HanXRQChr01  1166   1     
 6 HanXRQChr01  1335   0.352 
 7 HanXRQChr01  1404   1     
 8 HanXRQChr01  1429   0.352 
 9 HanXRQChr01  1450   1     
10 HanXRQChr01  1458   1     
# … with 6,740 more rows
```


**Side Question 1**:
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




