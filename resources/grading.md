#Course evaluation

The student evualtion in Biol 525D is based on the following:
* 50% participation
* 40% daily assignments
* 10% final shell script

###Participation
Attending class, participating in discussions and working through tutorials

###Daily assignments
As part of most lecture, there will be questions presented for you to work on and discuss in class. 
After the tutorial, we will select one question per topic and you will be responsible for answering that question. 
Answers should be less than one paragraph. 
Daily questions for topics 2-5 are due on _Monday August 22_, and for topics 6-10 on _Monday August 29_.
Please group your answers together into one email per week.
Answers should be sent to gregory.owens@alumni.ubc.ca

###Final shell script *Due Monday August 29th*
Students must produce a single shell script that does the following steps:
* Trims fastq files for quality.
* Aligns to a reference genome.
* Calls SNPs.
* Filters vcf output.
* Runs a principal component analysis.
* Outputs a pdf figure of the principal component analysis. 

Assume all programs and libraries are locally installed in the "/home/ubuntu/bin" directory. You can have additional scripts that are called by the shell script (i.e. an R script to run the PCA). 

#### R help
* To call an R script in command line use "R CMD BATCH Rscriptname.R"
* To print a pdf use the [pdf() command](http://www.cookbook-r.com/Graphs/Output_to_a_file/)
