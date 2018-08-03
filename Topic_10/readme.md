# Topic 10: Phylogenetics

The first step is to convert the VCF file into a fasta file. It is surprisingly hard to find a tool to do that, so I've supplied a script.

We only have variable sites, so we're going to use an ascertainment bias to correct for that. Typically phylogenetic programs assume you have all sites, including invariant ones.  
IQtree also has a [webserver](http://iqtree.cibiv.univie.ac.at/).

When using an ascertainment bias, every site needs to be variable. Annoyingly, when many phylogenetic programs treat heterozygous site as both homozygous types with some probability. Consequently, if an allele is only represented in its heterozygous state, there is some probability that the site is invariant, causing the program to crash. The script I provide filters out those sites. 

```

cd ~
zcat biol525d.snps.vcf.gz | perl /home/biol525d/bin/vcf2fasta_gaps.pl  > biol525d.snps.fasta



#Next we run it
iqtree -s biol525d.snps.fasta -st DNA -m TEST+ASC -nt 1

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
filename <- "Downloads/biol525d.snps.fasta.treefile"
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
#The samples are very clumpy because of the lack of variation. We don't have an outgroup species, so lets tree a midpoint root.

tree.midpoint <- midpoint(tree)
ggtree(tree.midpoint) +
  geom_tiplab() 
```

![](figure/ggtree-2.png)


```r
#How about we now label each population. I'm picking the node numbers off my phylogeny, and it may be different for your tree.
ggtree(tree.midpoint) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  geom_tiplab() +
  geom_cladelabel(node=111, label="Population 2", align=T, offset=.0001) +
  geom_cladelabel(node=199, label="Population 1", align=T, offset=.0001)

```

![](figure/ggtree-3.png)

```r
#Uh oh, labels are off the printed screen. Here's a work around
ggtree(tree.midpoint) +
  geom_text2(aes(subset=!isTip, label=node)) + 
  geom_tiplab() +
  geom_cladelabel(node=111, label="Population 2", align=T, offset=.01) +
  geom_cladelabel(node=199, label="Population 1", align=T, offset=.01) +
  geom_cladelabel(node=111, label="", align=T, offset=.03,color="white")

```

![](figure/ggtree-4.png)


```r
#We can  make up our labels
ggtree(tree.midpoint) +
  geom_tiplab() +
  geom_label2(aes(subset=(node==111), label='Robot Fish')) +
  geom_label2(aes(subset=(node==199), label='Robots')) +
  geom_cladelabel(node=111, label="", align=T, offset=.0001,color="white") + 
  geom_cladelabel(node=199, label="", align=T, offset=-.0002,color="white")

```

![](figure/ggtree-5.png)

### Coding challenge 1

Rename samples 050 and 051 to "turbo50" and "awesome51".

## Splitstree
Now we want to try out splitstree, a reticulate network phylogeny. Download it [here](http://ab.inf.uni-tuebingen.de/data/software/splitstree4/download/welcome.html). 
It should just run if you have the correct version of java. For newer macs you may have to install an older version of java, it should be a suggestion when you attempt to run the program.
Before we download the fasta file we have to replace all N with ?, because splitstree does not accept N.
```bash
cat biol525d.snps.fasta | sed s/N/?/g biol525d.snps.fasta > biol525d.snps.splitstree.fasta
```
Then transfer it to your laptop. Open Splitstree and select the biol525d.snps.splitstree.fasta file.
The first tree you see is a NeighbourNet, but you can also select many other algorithms or distance methods.


## Daily Questions:
1. What are two biological scenarios where you would want to use a reticulate tree, versus the standard non-reticulate tree, to represent the phylogeny?
2. What does it mean when something has 50% bootstrap support? What are two possible reasons that a node may have low support? Include one biological and one methodological reason.



