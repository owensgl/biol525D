# Topic 8: Population genomics

### Daily Questions:
1. For a site that is invariant in both populations (i.e. a locus with no variation), what is Fst?
2. Which two samples are the most extreme outliers (in either direction) for PC1?
3. What is the average Fst between the two populations in our example data? Hint, SNPrelate can calculate Fst.

The first thing we have to do is install some programs.

```bash
screen -r
cd bin
#ANGSD is a package for analyzing next gen sequencing data. It uses genotype uncertainty and can work off bam files directly.
wget http://popgen.dk/software/download/angsd/angsd0.912.tar.gz
tar xfz angsd0.912.tar.gz
cd htslib;make;cd ..
cd angsd
make HTSSRC=../htslib
cd ..
#NGSadmix
wget popgen.dk/software/NGSadmix/ngsadmix32.cpp 
g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

#FastStructure is similar to STRUCTURE but orders of magnitude faster. It has many dependencies

#Install anaconda, a python build for scientific computing
wget http://repo.continuum.io/archive/Anaconda2-4.1.1-Linux-x86_64.sh
bash Anaconda2-4.1.1-Linux-x86_64.sh
python=/home/ubuntu/anaconda2/bin/python
#Install GNU scientific library
wget http://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
tar -zxvf gsl-1.16.tar.gz
cd gsl-1.16
./configure
make
sudo make install

cd /home/ubuntu/bin
#Get faststructure
git clone https://github.com/rajanil/fastStructure

#Build the python extensions
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"
source ~/.bashrc
cd fastStructure/vars
$python setup.py build_ext --inplace
cd ..
$python setup.py build_ext --inplace
#Fast structure is now good to go! Remember to run it using the python version in anaconda
```
Now to convert between file types we're going to use PGDspider.

First download it to your home computer from [here](http://www.cmpg.unibe.ch/software/PGDSpider/#Download_and_Installation_Instructions)
We're going to create a spid file to help pgdspider convert a vcf to faststructure format. This [video](https://www.youtube.com/watch?v=I7hJvE0USxQ) demonstrates how.

We take the spid file created and transfer it to your server using cyberduck.

Next we download pgdspider to your server
```bash
java=/home/ubuntu/bin/jdk1.8.0_102/jre/bin/java
cd /home/ubuntu/bin/
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.0.3.zip
unzip PGDSpider_2.1.0.3.zip
#Now to convert the vcf file to faststructure format
cd ..
$java -jar bin/PGDSpider_2.1.0.3/PGDSpider2-cli.jar -inputfile biol525D.snps.vcf -inputformat VCF -outputfile biol525D.snps.str -outputformat STRUCTURE -spid VCF_to_structure.spid

#Now we run fast structure
mkdir faststructure
cd faststructure
for k in `seq 10`
do
$python /home/ubuntu/bin/fastStructure/structure.py -K $k --input=../biol525D.snps --output=biol525D --format=str
done
#This gives us a bunch of files for each K value, we'll use those to plot.
#The meanQ files give us admixture proportions which is what we care about.
#It's also useful to find out the best K value from the data and we can do that using faststructure
$python /home/ubuntu/bin/fastStructure/chooseK.py --input=biol525D
```
Now that we've run faststructure, its time to plot our results. First we download our faststructure folder to your laptop, then open Rstudio and continue onto the [next page](https://github.com/owensgl/biol525D/blob/master/Topic_8-9/plotting_structure.md)

