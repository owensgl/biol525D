# Calculating and plotting Fst

First we need to download and install vcftools
```bash
sudo apt-get install pkg-config
sudo apt-get install autoconf

cd /home/ubuntu/bin 
git clone https://github.com/vcftools/vcftools.git
cd vcftools/ 
./autogen.sh 
./configure	
make 
make install

#Next we use it to calculate Fst
/home/ubuntu/bin/vcftools/src/cpp/vcftools  \
--vcf biol525D.snps.vcf \
--weir-fst-pop samplelist.P1.txt --weir-fst-pop samplelist.P2.txt \
--out biol525D
#Note, samplelist.P1/2.txt are found on the github page. They are just lists of samples for each population.
#Take a look at the fst file, does it have reasonable data? Is it all missing data?
```
**Question 1**:
* Can you calculate overall Fst from this file and if so how would you do it using command line tools?

<details> 
  <summary>**Answer 1**  </summary>
   ```bash
    Fst is a ratio so calculating the overall values requires summing the numerator and denominator for each locus, which we don't have. 
    To calculate the mean value using command line you could use awk '{if ($3 != "-nan") sum += $3; n++ } END { print sum / n; }' biol525D.weir.fst
```
Next we'll take the Fst values and plot them in R. Transfer the biol525D.weir.fst file to your computer





