---
title: Course Evaluation
menuItem: Course Evaluation
layout: page
menuPosition: 2
---

The student evualtion in Biol 525D is based on the following:

* 50% participation
* 40% daily assignments
* 10% final shell script

### Participation
Attending class, participating in discussions and working through tutorials

### Daily assignments
As part of most lecture, there will be questions presented for you to work on and discuss in class. 
After the tutorial, we will select one question per topic and you will be responsible for answering that question. 
Answers should be less than one paragraph. 
Daily questions for all topics (2-10) are due on _Friday August 9_.
Please group your answers together into one email per week.
Answers should be sent to gregory.lawrence.owens@gmail.com

### Final shell script *Due Friday August 9th*
Students must produce a single shell script that does the following steps:
* Trims fastq files for quality.
* Aligns to a reference genome.
* Calls SNPs.
* Filters vcf output.
* Runs a principal component analysis.
* Outputs a pdf figure of the principal component analysis. 

Assume all programs and libraries are locally installed in the "/home/ubuntu/bin" directory. You can have additional scripts that are called by the shell script (i.e. an R script to run the PCA). The script should be able to run completely independently and not have any steps that require manual labor beyond the initial script call. 

### Grading
Scripts will be graded on the following criteria:
* Is the script generalizable? The script should work for any number of files with a specific naming scheme.
* Is the script annotated? The script should have comments throughout explaining what each step is doing.
* Is the script free of errors? The script should not have errors or skipped steps that would prevent it from working.

### R help
* To call an R script in command line use "R CMD BATCH Rscriptname.R"
* To print a pdf use the [ggsave() command](https://ggplot2.tidyverse.org/reference/ggsave.html)
