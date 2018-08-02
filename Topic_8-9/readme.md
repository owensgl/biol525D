# Topic 8: Population genomics

### Daily Questions:
1. For a site that is invariant in both populations (i.e. a locus with no variation), what is Fst?
2. What effect does missing data have on a PCA?
3. What is the average Fst between samples 1-50 and 51-100 in our example data? Hint, SNPrelate can calculate Fst.

A common first pass analysis is to use structure to look at clustering in your data. FastStructure is similar to STRUCTURE but orders of magnitude faster. We're going use that, but before that we have to convert our VCF to the bed format. We're going to use plink to do that.

```bash
cd ~/

#If you didn't make the biol525d.snps.vcf.gz file, you can download it from this github page. 
 plink --make-bed --vcf biol525d.snps.vcf.gz --out biol525d.snps --set-missing-var-ids @:# --double-id --allow-extra-chr

#Now we run fast structure with K values from 1 to 4.
#In this case, we're not explicitely filtering for linkage between SNPs (which you should), although subsetting to 10% of sites helps that somewhat.
mkdir faststructure
cd faststructure
for k in `seq 4`
do
python /home/biol525d/bin/fastStructure/structure.py -K $k --input=../biol525d.snps --output=biol525d
done
#This gives us a bunch of files for each K value, we'll use those to plot.
#The meanQ files give us admixture proportions which is what we care about.
#It's also useful to find out the best K value from the data and we can do that using faststructure
python /home/biol525d/bin/fastStructure/chooseK.py --input=biol525d
```
Now that we've run faststructure, its time to plot our results. First we download our faststructure folder to your laptop, then open Rstudio and continue onto the [next page](https://github.com/owensgl/biol525D/blob/master/Topic_8-9/plotting_structure.md)

