# Calculating and plotting Fst

First we need to download and install vcftools. This is back on your ubuntu server.
```bash
sudo apt-get install pkg-config
sudo apt-get install autoconf

cd ~/bin 
git clone https://github.com/vcftools/vcftools.git
cd vcftools/ 
./autogen.sh 
./configure	
make 
sudo make install

#Next we use it to calculate Fst
#You will have to download from the github page samplelist.P1.txt and samplelist.P2.txt
cd ~
vcftools  \
--vcf Biol525D.snps.vcf \
--weir-fst-pop samplelist.P1.txt --weir-fst-pop samplelist.P2.txt \
--out Biol525D
#Note, samplelist.P1/2.txt are found on the github page. They are just lists of samples for each population.
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
##       groupIV      groupVII        groupI       groupIX       groupXX 
##          1234          1201          1081           806           806 
##       groupII     groupVIII       groupXI      groupXII     groupXIII 
##           793           748           652           651           645 
##     groupXVII    groupXVIII       groupVI      groupIII      groupXVI 
##           600           594           591           554           554 
##      groupXIV        groupX      groupXIX       groupXV      groupXXI 
##           546           542           532           482           482 
##        groupV   scaffold_27   scaffold_37   scaffold_99   scaffold_74 
##           428           144            71            55            45 
##   scaffold_67   scaffold_47   scaffold_48   scaffold_68  scaffold_122 
##            44            37            37            35            29 
##   scaffold_54  scaffold_168  scaffold_112  scaffold_115  scaffold_184 
##            28            20            19            19            19 
##  scaffold_151   scaffold_95  scaffold_114  scaffold_126  scaffold_149 
##            16            16            14            14            13 
##   scaffold_56  scaffold_120  scaffold_165  scaffold_214   scaffold_84 
##            13            12            12            12            12 
##  scaffold_135   scaffold_61  scaffold_111  scaffold_157  scaffold_169 
##            10            10             9             9             9 
##  scaffold_211  scaffold_256   scaffold_69   scaffold_89   scaffold_90 
##             9             9             9             9             8 
##   scaffold_98  scaffold_180  scaffold_202  scaffold_213  scaffold_300 
##             8             7             7             7             7 
##  scaffold_354  scaffold_101  scaffold_132  scaffold_139  scaffold_150 
##             7             6             6             6             6 
##  scaffold_156  scaffold_159  scaffold_177  scaffold_181  scaffold_197 
##             6             6             6             6             6 
##  scaffold_206  scaffold_363   scaffold_58   scaffold_88  scaffold_129 
##             6             6             6             6             5 
##  scaffold_137  scaffold_161  scaffold_163  scaffold_175  scaffold_216 
##             5             5             5             5             5 
##  scaffold_238  scaffold_246  scaffold_270  scaffold_273  scaffold_281 
##             5             5             5             5             5 
##  scaffold_545  scaffold_740   scaffold_80 scaffold_1075  scaffold_121 
##             5             5             5             4             4 
##  scaffold_128  scaffold_133  scaffold_152  scaffold_182  scaffold_204 
##             4             4             4             4             4 
##  scaffold_208  scaffold_225  scaffold_249  scaffold_369       (Other) 
##             4             4             4             4           121
```

```r
#This strips the "group" from the chromosome name
data$CHROM <- gsub("group", "", data$CHROM)

#This converts it to numeric
data$CHROM <- as.numeric(as.roman(data$CHROM))
#This removes scaffolds
data %>% filter(!is.na(CHROM)) -> data

#Lets also remove values that are NA
data %>% filter(WEIR_AND_COCKERHAM_FST != "NaN") -> data

#It's important to make sure that the chromosomes are numeric instead of character
data$CHROM <- as.numeric(data$CHROM)

#Lets do a basic plot using the manhattan tool in qqman. This is generally designed for plotting pvalues from GWAS, but it works here.
manhattan(data, chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST", snp="POS",
          logp=FALSE,ylab="Fst")
```

![](figure/fst1-1.png)


**Question 2**:
* Why are there so few possible Fst values? Why are there so few that have zero Fst?

**Coding challenge**:
* Do a manhattan plot in ggplot or base R. As a bonus, make it look nice.




