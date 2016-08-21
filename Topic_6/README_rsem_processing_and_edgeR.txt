

####################################################################################
###### RNAseq analysis exercise PART I: RSEM Processing
####################################################################################

Steps for running RSEM:

see documentation here for more info:

http://deweylab.biostat.wisc.edu/rsem/

Typing in "rsem-calculate-expression" or any of the other commands without any arguments will bring up a help screen. In all the RSEM commands below, I have put the syntax for the command after "##command", including terms that look like "<something>". When you see this, the "something" will tell you what should go in that argument, and the "<>" characters are there to indicate to you that it's an argument, and should not be included when you actually run it. There is also an example command, indicated by "##example"


Step 1. Install Bowtie.

##command:
sudo apt-get install bowtie2


Step 2. Install RSEM (~/programs/RSEM-1.2.31.tar.gz).

First unpack.

##commands:
gunzip -d RSEM-1.2.31.tar.gz

tar -xvf RSEM-1.2.31.tar

Then compile the program. 

##commands:
cd RSEM-1.2.31
sudo make
sudo make install


RSEM also needs Rscript to run:

##command:
sudo apt-get install r-base-core


Run the remaining commands from the ~/data/Topic6/ directory. 

Step 3. Build the transcript-to-gene map using RSEM. This associates each gene isoform to a gene name.

##command:
extract-transcript-to-gene-map-from-trinity <fasta_reference> <output_map_name>

##example:
extract-transcript-to-gene-map-from-trinity Pine_reference_rnaseq_reduced.fa Pine_reference_rnaseq_reduced_map


Step 2. Prepare the reference, so that RSEM can use it:

##command:
rsem-prepare-reference --transcript-to-gene-map <map_name_from_step1> <fasta_reference> <outputname>

##example:
rsem-prepare-reference --transcript-to-gene-map Pine_reference_rnaseq_reduced_map Pine_reference_rnaseq_reduced.fa Pine_reference_rnaseq_reduced_ref


Step 3. Prepare a bowtie2 index of the fasta file with the same output_name as the reference output in step #2. Bowtie is the alignment program that will align the reads to the reference:

##command:
bowtie2-build -f <fasta_reference> <output_name>

##example:
bowtie2-build -f Pine_reference_rnaseq_reduced.fa Pine_reference_rnaseq_reduced_ref


Step 4. Calculate expression. Previous versions of RSEM specified the "--no-polyA" flag to be used if the data have already been cleaned to remove polyA tails, but this is now the default, so poly-A tails are not added unless you specify. If you run RSEM in the same directory that the fq files (~/data/Topic6) are in, and do not provide an output path, the output will appear in that directory.


##command:
rsem-calculate-expression --bowtie2 --paired-end <fastq_R1> <fastq_R2> <name_of_prepared_ref_step1> <output_name>

##example:
rsem-calculate-expression --bowtie2 --paired-end PmdT_147_100k_R1.fq PmdT_147_100k_R2.fq Pine_reference_rnaseq_reduced_ref PmdT_147_rsem


Examine the output of RSEM, both the output to the screen and the files that it created. You should be able to see the percentage of reads that were aligned correctly, among other summary statistics. Refer to the online manual for a description of the files. The description of the  see RSEM-1.2.31/cnt_file_description.txt

#################
Recap: 

Use "ls" to discover what the output is called, and use "less" or "head" to see what the files look like.

################

Now rerun this command on the other two paired .fq libraries, generating three RSEM alignments. 

Step 5. To plot some quality statistics, run the following command, replacing the <NAME> with the appropriate filename prefix, without the .stat ending (and removing the <> characters):

##command:
rsem-plot-model <name_of_aligned_output_from_step4> outputname.pdf

##example:
rsem-plot-model PmdT_147_rsem PmdT_147_rsem.pdf

From the RSEM manual:
The plots generated depends on read type and user configuration. It may include fragment length distribution, mate length distribution, read start position distribution (RSPD), quality score vs observed quality given a reference base, position vs percentage of sequencing error given a reference base and histogram of reads with different number of alignments.

fragment length distribution and mate length distribution: x-axis is fragment/mate length, y axis is the probability of generating a fragment/mate with the associated length

RSPD: Read Start Position Distribution. x-axis is bin number, y-axis is the probability of each bin. RSPD can be used as an indicator of 3’ bias

Quality score vs. observed quality given a reference base: x-axis is Phred quality scores associated with data, y-axis is the “observed quality”, Phred quality scores learned by RSEM from the data. Q = –10log_10(P), where Q is Phred quality score and P is the probability of sequencing error for a particular base

Position vs. percentage sequencing error given a reference base: x-axis is position and y-axis is percentage sequencing error

Histogram of reads with different number of alignments: x-axis is the number of alignments a read has and y-axis is the number of such reads. The inf in x-axis means number of reads filtered due to too many alignments

You can download these to your computer to take a look at the graphs.

Step 6. To make a table out of the individual library expression files that you have created, use the custom perl script called "add_RSEM_data_to_table.pl". To run this, you will first have to create a list of the contig names that were in your reference, along with a list of the .genes.results files created by RSEM that you want to put together. The list of names in the reference will have to match the names that you see in the RSEM output files that end in ".genes.results". See if you can figure out how to do these operations using simple bash commands. The commands to do so are listed at the bottom of this document if you run into trouble. Then run the perl script:


##command:
perl add_RSEM_data_to_table.pl <list_to_add> <gene_names> <suffix_for_output>

##example:
perl ~/scripts/add_RSEM_data_to_table.pl list_to_add.txt gene_names.txt _expression_table.txt

You can see that this outputs a table that is readable in R, with one row for each gene and one column for each individual, with the expected counts from the 5th column printed out. If you wish to modify what is printed to the file, it can be done by editing line 81.

Question 1. What is the expected count of comp996_c0 for each individual?

Question 2. What expression measure would you use to compare gene expression between different genes and why? Is it appropriate to compare the raw expression counts? Can you get more appropriate data from RSEM? Discuss this with a group of four and be prepared to share your ideas with the class.


####################################################################################
###### RNAseq analysis exercise PART II: EdgeR analysis
####################################################################################

There are two data files:
cold_hot_expression.txt
cold_hot_mwh_expression.txt

They have been created by using RSEM to align libraries to a lodgepole pine reference and represent a subset of the data published in Yeaman et al. (2014; New Phytologist). Individuals were grown in a common garden and then exposed to conditions that were either hot (H), cool (C), or mild and wet with a heat treatment (MWH; see paper for details). The reference dataset has been trimmed so that it is smaller and easier to work with for this workshop.

The script "process_hot_cold_expression.R" will show you some of the basic steps for analyzing expression data. For more details, the EdgeR user guide (included in the folder here) provides an excellent resource with well worked examples. Use these examples to play around with the data and try to answer the following questions:

Question 1. How many genes are differentially expressed by treatment in the simple contrast of C vs H (using dataset "cold_hot_expression.txt")? How does the choice of FDR cutoff or p-value affect this number? 

Question 2. How many genes are differentially expressed in the three-way contrast (using "cold_hot_mwh_expression.txt")? Which treatment is driving differential expression here? How do you know?

Question 3. How much does model fitting with common dispersion vs. tagwise dispersion affect the answers you get from the data? (think in terms of the number of DE genes, the evidence for a single gene, etc.)


##### STEPS for bash commands to prepare input files for add_RSEM_data_to_table.pl 

List all files with .genes.results in their name (you may want to delete some if you've made more copies than you should have during testing, or chain multiple "| grep" commands together)

ls | grep genes.results > list_to_add.txt


Find all lines that have ">" in them, which are the contig names. Then pass these to sed and strip off the ">" character, and save it to a file:

grep ">" Pine_reference_rnaseq_reduced.fa | sed 's/>//' > gene_names.txt






 