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
sudo apt-get install python-setuptools
sudo apt-get install python-pip
git clone http://github.com/numpy/numpy.git numpy
cd numpy
sudo python setup.py install
cd ..
sudo rm -R numpy
