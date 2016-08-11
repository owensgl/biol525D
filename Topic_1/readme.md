# Topic 1: Setting up programs

In this course we need several programs. 
Programs to install: 
* [Rstudio](https://www.rstudio.com/products/rstudio/download2/)
* [github desktop](https://desktop.github.com/)
* [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) *For windows
* [Cyberduck](https://cyberduck.io/?l=en)

Log into AWS server using ssh (terminal or putty)

```bash
#Update the OS
sudo apt-get update

#Install essential programs
sudo apt-get -y install cmake
sudo apt-get -y install build-essential
sudo apt-get -y install parallel
sudo apt-get -y install zip
sudo apt-get -y install zlib1g-dev
sudo apt-get -y install libncurses5-dev
sudo apt-get -y install git
cd bin
 wget --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u102-b14/jdk-8u102-linux-x64.tar.gz
tar xzf jdk-8u102-linux-x64.tar.gz



