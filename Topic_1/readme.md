---
title: "Topic 1: Logging into the server."
permalink: /Topic_1/
topickey: 1
topictitle: Servers
---

We are working on servers provided by Compute Canada, hosted on the [WestCloud](https://www.computecanada.ca/research-portal/national-services/compute-canada-cloud/) system. The only task now is to install a few programs we'll need later on, and login to your account.

Accompanying material:
* [Slides](./Topic 1.pdf)

Programs to install:
* [Rstudio](https://www.rstudio.com/products/rstudio/download2/)
* [Cyberduck](https://cyberduck.io/?l=en)
* [iTerm2](https://www.iterm2.com/)

Log into WestCloud server using ssh (terminal). You will find your username and ip address in this [table](https://docs.google.com/spreadsheets/d/1v7k2-XtfiwOoQ3iZHnJyqVXsxgekVGXEtnFIvdk7aqU/edit?usp=sharing).
```bash

ssh  -v -i ~/.ssh/id_rsa <your_user_name>@<your_server_ip>
```

Your screen should look like this:

![](terminal.jpeg "Terminal"){:width="100%"}


