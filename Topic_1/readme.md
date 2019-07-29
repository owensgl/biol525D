---
title: "Topic 1: Logging into the server."
permalink: /Topic_1/
topickey: 1
topictitle: Servers
---

We are working on servers provided by Compute Canada, hosted on the
[WestCloud](https://www.computecanada.ca/research-portal/national-services/compute-canada-cloud/)
system. The only task now is to install a few programs we'll need
later on, and login to your account.

**Accounting**: We have created an account for each student. You may
find your account information in this
[table](https://docs.google.com/spreadsheets/d/1v7k2-XtfiwOoQ3iZHnJyqVXsxgekVGXEtnFIvdk7aqU/edit?usp=sharing). Take
note of your username, and server ip address. In the instructions, if
you see the placeholders `serveruser` and `serverhost` replace them
with your own information. You should have received your server password separately.


Accompanying material:
---------------------

* [Slides](./Topic 1.pdf)



Software: _(can be completed early)_{: style="color: darkred"}
-------------------

Your local workstation (not the server) will need to be fitted with some
software for programming in R, and connecting to servers. The servers
will be fitted with the rest of the software already.

On your computer, you will need:

1. [Rstudio](https://www.rstudio.com/products/rstudio/download2/), the de-facto Integrated Development Environment (IDE) for R.
   The free edition will suffice for this course.

1. Tool to transfer files over SFTP (secure file transfer protocol) or
   SCP (secure copy). We can use these tools to browse and transfer files
   back and forth between your workstation and the server. Get *one* of the following options:

   * _Graphical option 1_: [Cyberduck](https://cyberduck.io/?l=en) (Mac &
    Windows). This is a graphical tool that will allow you to transfer
    files to/from a remote storage location (+dragndrop), and edit remote
    files with a local editor. We will also cover how to perform these
    operations from the command line using other readily-availble
    tools.

   * _Graphical option 2_:
   [MobaXTerm](https://mobaxterm.mobatek.net/download-home-edition.html)
   (Windows only). MobaXTerm combines a graphical terminal emulator, a
   graphical file-transfer tool, an X11 server (to use graphical
   applications from the remote-side), and more. Comes in two formats:
   "portable", and "installer/msi". The difference: the portable is
   installed simply by uncompressing the contents of the zip in a
   folder of your choice (this is useful if you want to carry it
   around on a USB stick). The "installer edition" will run an installation
   wizard and get added to your start menu and windows environment.

   * _Command line option_1_: pscp.exe/psftp.exe from the [PuTTy](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) suite (Windows only). Install the MSI version to get all the utilities at once. It comes in 32bit and 64bit editions. Unless you know your workstation has a 32bit cpu or operating system (which is uncommon in 2019), pick 64bit.

   * _Command line option 2_: Other platforms (Mac & GNU/Linux) have
     sftp and scp tools installed out-of-the-box. We will also teach
     you how to transfer files from the command-line using these. This
     knowledge will come in handy for scripting.

1. A (good) terminal emulator. To interact with command line programs, you will need to type a terminal emulator, and parse results.

   If you don't have a terminal emulator, we recommend you install one
   from the list that follows (there are many more out there):

    - [iTerm2](https://www.iterm2.com/) (Mac Only). Fuller-featured replacement for the default Terminal application (see feature list). The built-in terminal application which comes with MacOs is also decent.

    - [MobaXTerm](https://mobaxterm.mobatek.net/download-home-edition.html) (Windows). This is a one-stop shop for working with remote computers from windows. It has a good terminal, and it also doubles as a file transfer tool.

    - putty.exe from the [PuTTy](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html) suite (Windows). This is a second choice if Moba doesn't work. Putty is a "classic" terminal emulator on Windows, because it is stable, and simple to use. Get the .msi (installer) version which has all of the utilities. It comes in 32bit and 64bit editions. Unless you are working with a 32bit cpu or operating system (which is uncommon in 2019), pick 64bit.

    - For GNU/Linux environments (e.g. Ubuntu). Your desktop environment's default terminal application will work just fine, e.g. gnome-terminal, or xfce-terminal.

   Features of a good terminal emulator (things to look for):
     1. Options for scrolling and searching.
     1. Reliable text selection, copy-paste, both in, and out.
     1. Sensible translation of keypresses (especially true of <kbd>ALT</kbd>, <kbd>CTRL</kbd>, <kbd>⌘cmd</kbd>, <kbd>PGUP</kbd>, <kbd>PGDN</kbd>, <kbd>DEL</kbd>, <kbd>BCKSPC</kbd>, <kbd>⊞ Win</kbd>, <kbd>FN</kbd>, and arrow keys, whose meaning can be open to different interpretation across different platforms).
     1. Proper display of characters from different unicode planes (natural languages, emojis, etc.)
     1. Configurable appearance (font sizes, colors, and transparency).
     1. Multiple tabs/sessions in the same window.

   > *Note*: I would advise *against* using Window's CMD prompt for
   this purpose because it's too minimal. In a pinch, it can do, but
   it's worthwhile to spend a few minutes on getting a better
   emulator with better defaults.
   >
   *Note 2*: If you're curious why it's called an "emulator", you may read about their ancestor: ["real" terminals](https://en.wikipedia.org/wiki/Computer_terminal#Text_terminals).


Configuration: _(can be completed early)_{: style="color: darkred"}
---------------

This section assumes that you've already followed the steps to install
software above for your environment. We have provided steps to
configure rapid authentication between your computer and the
server. The steps are meant to be done in sequence and should take
about 15 minutes to complete once the software above is installed.

Configuration steps, in sequence:

1. Generate a key for authenticating to your assigned server. Follow steps in [Generate A Key](./generate_a_key).

2. Load the key in your ssh-agent. Details in [Configure-SSH-Agent](./configure_ssh_agent)

3. Configure the server account to authenticate with your new key. Take note of your assigned server hostname (IP), username, and password and follow details in [Finalize Tool Config](./finalize_tool_config). This will let you login to the server without remembering usernames, ips, or keys.


Testing it all: _(can be completed early)_{: style="color: darkred"}
---------------

Log into WestCloud server using ssh (terminal). Recall that you will find your username and ip address in this [table](https://docs.google.com/spreadsheets/d/1v7k2-XtfiwOoQ3iZHnJyqVXsxgekVGXEtnFIvdk7aqU/edit?usp=sharing). We use `serveruser` and `serverhost` as placeholders for this info.

The information presented by the servers when you connect to them is shown here in [Fingerprints](./fingerprints).

* For Mac, GNU Linux users:

  ```bash

  ssh -v <serveruser>@<serverhost>

  # if you have followed the key authentication steps, you need only to specify
  # the host alias you chose from your `~/.ssh/config`, e.g. b525:

  ssh -v b525
  ```
* MobaXTerm users:

  If you have configured a new user session, as specified in the
  key authentication pages, you can simply find your session
  on the left pane of Moba (the vertical-text tab "Sessions"), and
  double click your chosen session.

  If not, you can open a moba local terminal and do:

  ```
  # from the local terminal
  ssh -v <serveruser>@<serverhost>

  # or just this if you have configured ~/.ssh/config
  ssh -v b525
  ```

* Putty users:

  If you've followed the key authentication steps, you should have a user session you can open
  directy from the Pageant right-click menu.

  Otherwise, create a new putty session:

  - In the menu, use your provided account info for Host Name (your `serverhost`), Connection->Data->Login username (your `serveruser`).
  - If you have followed the key authentication steps, you will not have to enter
    passwords, and will be logged in immediately. Otherwise, provide your password at the prompt (your `serverpass`).

Your screen should give you a prompt on the server, and should print some server text similar to this (ignoring visual terminal differences aside):

![](terminal.jpeg "Terminal"){:width="100%"}


