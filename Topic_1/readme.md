# Topic 1: Setting up programs

In this course we need several programs. 
Programs to install: 
* [Rstudio](https://www.rstudio.com/products/rstudio/download2/)
* [MobaXterm](http://mobaxterm.mobatek.net/download.html) *For windows
* [Cyberduck](https://cyberduck.io/?l=en)

Log into WestCloud server using ssh (terminal or MobaXterm). You will receive an IP address through email. 
```bash
ssh  -v  ubuntu@<your_server_ip>
```


Now that you've logged into the server, next we have to install and update some programs.

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
mkdir bin
cd bin
 wget --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u102-b14/jdk-8u102-linux-x64.tar.gz
tar xzf jdk-8u102-linux-x64.tar.gz
```


