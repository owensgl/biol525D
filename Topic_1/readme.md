# Topic 1: Setting up programs

In this course we need several programs. 
Programs to install: 
* [Rstudio](https://www.rstudio.com/products/rstudio/download2/)
* [github desktop](https://desktop.github.com/)
* [Putty](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html) *For windows
* [Cyberduck](https://cyberduck.io/?l=en)

Log into AWS server using ssh (terminal or putty).
#### Mac OS or Linux

* Download Server_access.zip to Downloads folder
* Unzip Server_access.zip (automatic on mac)
* Open test_server_ssh.bash and replace <INSERT IP HERE> with IP from your server on [list](https://docs.google.com/spreadsheets/d/1k3o-g60c_3Parf0HDOfSnkjwemwYatn3UJLrNISMUNs/edit?usp=sharing)
* open terminal
```bash
cd Downloads
cd Server_access
bash ./test_server_ssh.bash
```
#### Windows
* Download putty and puttygen from [here](http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html)
* Download Server_access.zip to Downloads folder
* Unzip Server_access.zip
* Follow [this tutorial](https://support.rackspace.com/how-to/logging-in-with-an-ssh-private-key-on-windows/) to log into the server. The private key is biol525D.pem from Server_access folder, your IP address is from this [list](https://docs.google.com/spreadsheets/d/1k3o-g60c_3Parf0HDOfSnkjwemwYatn3UJLrNISMUNs/edit?usp=sharing).


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
cd bin
 wget --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" http://download.oracle.com/otn-pub/java/jdk/8u102-b14/jdk-8u102-linux-x64.tar.gz
tar xzf jdk-8u102-linux-x64.tar.gz



