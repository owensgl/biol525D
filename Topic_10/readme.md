# Topic 10: Phylogenetics

The first step is to convert the formatted.tab file into a fasta file. The script will convert the tab style file that we created to a fasta, for heterozygous sites it picks one random allele. IQtree and Splitstree both have seem to ignore ambiguous bases, so we're not going to use them.

We also have to remove invariant sites, so that we can use an ascertainment bias in IQtree. 
IQtree also has a [webserver](http://iqtree.cibiv.univie.ac.at/).

```
cd ~
#Note: Biol525D.snps.formatted.tab comes from topic 7
cat Biol525D.snps.formatted.tab | perl /mnt/data/bin/SNPtable2fasta_pickoneallele.pl | perl /mnt/data/bin/fasta2removeinvariant.pl >Biol525D.snps.fasta


#Next we should install IQ-tree
cd ~/bin
wget https://github.com/Cibiv/IQ-TREE/releases/download/v1.5.5/iqtree-omp-1.5.5-Linux.tar.gz
tar -xzf iqtree-omp-1.5.5-Linux.tar.gz

#Next we run it
bin/iqtree-omp-1.5.5-Linux/bin/iqtree-omp -s Biol525D.snps.fasta -st DNA -m TEST+ASC

#This produces several output files, including a log and a couple different versions of the treefile. The -m TEST command does a model test and selects the best substition model. 
```
In the next step, we're going to use R to visualize our tree using ggtree. To do that you need to download "Biol525D.snps.fasta.treefile" to your laptop. Another way to visualize a tree is [Figtree](http://tree.bio.ed.ac.uk/software/figtree/).


```r
#First we install some packages
source("http://bioconductor.org/biocLite.R")
biocLite("ggtree")
install.packages("phytools")
```

```r
#Then load some libraries
library(ggtree)
library(phytools)
library(phangorn)
```



```r
#Then load in our data
filename <- "/Downloads/Biol525D.snps.fasta.treefile"
tree <- read.tree(filename)

#Lets take a look
ggtree(tree, layout="unrooted") +
  #This labels nodes by their number, so we can work with them.
  geom_text2(aes(subset=!isTip, label=node)) + 
  #This labels tips.
  geom_tiplab() 
```

![](figure/ggtree-1.png)


```r
#Since we have three samples, we don't have a known root. 
#In this case we should do a midpoint root between the two populations

tree.midpoint <- midpoint(tree)
ggtree(tree.midpoint) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  geom_tiplab() 
```

![](figure/ggtree-2.png)


```r
#How about we now label a group
ggtree(tree.midpoint) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  geom_tiplab() +
  geom_cladelabel(node=5, label="Population 1", align=T, offset=.0001)
```

![](figure/ggtree-3.png)

```r
#Uh oh, labels are off the printed screen. Here's a work around
ggtree(tree.midpoint) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  geom_tiplab() +
  geom_cladelabel(node=5, label="Population 1", align=T, offset=.0001) +
  geom_cladelabel(node=5, label="", align=T, offset=.0002,color="white")
```

![](figure/ggtree-4.png)


```r
#We can  make up our labels
ggtree(tree.midpoint) +
  geom_tiplab() +
  geom_label2(aes(subset=4, label='Robot Fish')) +
  geom_label2(aes(subset=5, label='Robots')) +
    geom_cladelabel(node=5, label="", align=T, offset=.0001,color="white") + 
    geom_cladelabel(node=5, label="", align=T, offset=-.0002,color="white")
```

![](figure/ggtree-5.png)

### Coding challenge 1

Replace the names of the three samples with the names of the three stooges. 

## Splitstree
Now we want to try out splitstree, a reticulate network phylogeny. Download it [here](http://ab.inf.uni-tuebingen.de/data/software/splitstree4/download/welcome.html). 
It should just run if you have the correct version of java. For newer macs you may have to install an older version of java, it should be a suggestion when you attempt to run the program.
Before we download the fasta file we have to replace all N with ?, because splitstree does not accept N.
```bash
cat Biol525D.snps.fasta | sed s/N/?/g Biol525D.snps.fasta > Biol525D.snps.splitstree.fasta
```
Then transfer it to your laptop. Open Splitstree and select the biol525D.snps.splitstree.fasta file.
The first tree you see is a NeighbourNet, but you can also select many other algorithms or distance methods.

The first tree was very boring, so lets try a bigger dataset. Download example.fa from the github page and then run that in Splitstree. Use the zoom controls to make a readable tree. Export the figure to svg for further editting. With a tree like this, how would you format it to make it understandable in a paper?

## Daily Questions:
1. What are two biological scenarios where you would want to use a reticulate tree, versus the standard non-reticulate tree, to represent the phylogeny?
2. What does it mean when something has 50% bootstrap support? What are two possible reasons that a node may have low support? Include one biological and one methodological reason.



