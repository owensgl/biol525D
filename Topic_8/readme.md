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
sudo apt-get install python-dev
#Cython...
wget https://pypi.python.org/packages/c6/fe/97319581905de40f1be7015a0ea1bd336a756f6249914b148a17eefa75dc/Cython-0.24.1.tar.gz#md5=890b494a12951f1d6228c416a5789554
tar xzf Cython-0.24.1.tar.gz
cd Cython-0.24.1/
sudo python setup.py install
cd ..
#Numpy...
git clone http://github.com/numpy/numpy.git numpy
cd numpy
sudo python setup.py install
cd ..
#Scipy...
git clone http://github.com/scipy/scipy.git scipy
cd scipy
sudo python setup.py install
cd ..
#Now clean up installers
sudo rm -R numpy
sudo rm -R scipy
sudo rm -R 

#Install anaconda 
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



