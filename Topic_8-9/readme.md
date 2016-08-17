# Topic 8: Population genomics

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

#Install GNU scientific library
wget http://gnu.mirror.vexxhost.com/gsl/gsl-latest.tar.gz
tar -zxvf gsl-latest.tar.gz
cd gsl-1.16
./configure
make
sudo make install

#Get faststructure
git clone https://github.com/rajanil/fastStructure

#Build the python extensions
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
export CFLAGS="-I/usr/local/include"
export LDFLAGS="-L/usr/local/lib"
source ~/.bashrc
cd fastStructure/vars
/home/ubuntu/anaconda2/bin/python setup.py build_ext --inplace
cd ..
/home/ubuntu/anaconda2/bin/python setup.py build_ext --inplace
#Fast structure is now good to go! Remember to run it using the python version in anaconda

#Now to convert between file types we're going to use PGDspider
#First download it to your home computer from [here](http://www.cmpg.unibe.ch/software/PGDSpider/#Download_and_Installation_Instructions)
#We're going to create a spid file to help pgdspider convert a vcf to faststructure format. This [video](https://www.youtube.com/watch?v=I7hJvE0USxQ) demonstrates how.
#We take the spid file created and transfer it to your server using cyberduck.
#Next we download pgdspider to your server
cd /home/ubuntu/bin/
wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.1.0.3.zip
unzip PGDSpider_2.1.0.3.zip



