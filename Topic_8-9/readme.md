---
title: "Topics 8: Population genomics"
permalink: /Topic_8-9/
topickey: 8
topictitle: "Population genomics"
---

## Accompanying material
* [Slides](./Topic 8.pdf)

## Daily Questions:
1. For a site that is invariant in both populations (i.e. a locus with no variation), what is Fst?
2. What effect does missing data have on a PCA?
3. What is the average Fst between ARG and ANN samples in our dataset? Hint, SNPrelate can calculate Fst.

Last topic we called variants across the three chromosomes. If you look at the VCF, you'll notice there are a lot of sites only genotyped in a small subset of the samples. This can happen with lower overall read depth (in this case this is whole genome sequencing at ~7X depth), but due to other factors like divergence been sample and reference. We also have indels, and SNPs with more than two alleles. Many programs strictly require biallelic sites so lets first filter the VCF to a smaller set of usable sites.
We're going to use _bcftools_ a very fast program for processing and filtering VCF files. Here are what we want to filter for:
* At least 8/10 samples genotyped (which are 16 alleles since they're diploid.
* Only one alternate allele.
* No indels.
* At least 2 copies of the alternate allele

```bash
bcftools  -c 2 -i 'INFO/AN >= 16' -m2 -M2 -v snps vcf/full_genome.vcf.gz -O z > vcf/full_genome.filtered.vcf.gz
```


A common first pass analysis is to use structure to look at clustering in your data. FastStructure is similar to STRUCTURE but orders of magnitude faster. We're going use that, but before that we have to convert our VCF to the bed format. We're going to use plink to do that.


```bash
cd ~/

    # See note above.
    /mnt/bin/plink-1.07-x86_64/plink --make-bed \
          --noweb \
    	  --vcf vcf/full_genome.filtered.vcf.gz \
          --out vcf/full_genome.filtered \
	  --set-missing-var-ids @:# \
	  --double-id \
	  --allow-extra-chr

# Now we run fast structure with K values from 1 to 4.
# In this case, we're not explicitely filtering for linkage between SNPs (you should),
# although subsetting to 10% of sites helps that somewhat.
mkdir faststructure
cd faststructure
for k in `seq 4`; do
  python /home/biol525d/bin/fastStructure/structure.py -K $k \
         --input=../biol525d.snps \
	 --output=biol525d
done

# This gives us a bunch of files for each K value, we'll use those to plot.
# The meanQ files give us admixture proportions which is what we care about.
# It's also useful to find out the best K value from the data and we can do
# that using faststructure
python /home/biol525d/bin/fastStructure/chooseK.py --input=biol525d
```

Now that we've run faststructure, its time to plot our results.
First we download our faststructure folder to your laptop, then open Rstudio
and continue onto the [next page](./plotting_structure.md).

