---
title: "Topic 9: Plotting FastStructure results in R"
topickey: 9
topictitle: "Plotting in R"
datafiles: [1,2,3,4]
---

Make sure you have downloaded your faststructure directory to your Downloads folder using cyberduck or scp.

``` r
#First we install some packages

install.packages("tidyverse")
install.packages("reshape")


#Then we load those packages.

library(reshape)
library(tidyverse)
```


Save these files to your Downloads (these files are also in the course repository under ./Topic_8-9/):
  - Sample info [biol525d.popinfo.txt](biol525d.popinfo.txt)
  - Data files:
    {% for i in page.datafiles %}
    - [biol525d.{{i}}.meanQ](./biol525d.{{i}}.meanQ), [biol525d.{{i}}.meanP](./biol525d.{{i}}.meanP)
    {% endfor %}

``` r
#"biol525d.popinfo.txt" has information on our samples.
sampleinfo <- read_tsv("Downloads/biol525d.popinfo.txt")
#data name and directory
name <- "Downloads/faststructure/biol525d"
#We're going to loop through each k value, so we need a dataframe to save those values
data.full <- data.frame(name=character(),
                        population=character(),
                        lat=numeric(),
                        long=numeric(),
                        variable=character(),
                        value=numeric(),
                        k=numeric())
#Now we loop through each K output
for (k in 1:4){
  #Load the data file
  data <- read.table(paste(name, k, "meanQ",sep="."),colClasses="numeric")
  #Label the columns (one for each group)
  Qnames <-paste("Q",seq(1:k),sep = "")
  colnames(data) <- Qnames
  #Add sample info to Q values
  data.info <- cbind(data, sampleinfo)
  #Melt the data into long format
  data.melt <- melt(data.info, id.var=c("name","population","lat","long"))
  #We have to make sure to include the K value for each file
  data.melt$k <- k
  #Now rbind it to the data frame outside the loop
  data.full <- rbind(data.full, data.melt)
}
#Lets try plotting for k=2
data.full %>% filter(k == 2) %>% #This selects only the k=2 lines out of the full set
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack")
```

![](figure/structure1-1.png)

``` r
#From this, we can easily scale up to all k values using facet
data.full %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=2)
```

![](figure/structure1-2.png)

``` r
#How about if we want to order samples by the reverse of latitude.
data.full %>%
  mutate(name = factor(name, levels = name[order(-lat)])) %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=2)
```


![](figure/structure1-3.png)

``` r
#The order works, but lets try to make it look nicer

data.full %>%
  mutate(name = factor(name, levels = name[order(-lat)])) %>%
  ggplot(.,aes(x=name,y=value,fill=factor(variable))) +
  geom_bar(stat = "identity",position="stack") +
  facet_wrap(~k, nrow=5) +
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        strip.background = element_blank()) +
  ylab("q-value")+
  guides(fill = guide_legend(title = "Group", title.position = "left"))
```


![](figure/structure1-4.png)

Plotting challenge 1
--------------------

-   For K = 2, plot the faststructure results and divide your samples by populations. Furthermore, make group 1 red and group 2 blue. Title the plot "*Helianthus* is great!"



Now lets move onto [Principal Component Analysis](./pca.md)
