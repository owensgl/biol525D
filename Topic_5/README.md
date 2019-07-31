---
title: "Topic 5: Assembly"
permalink: /Topic_5/
topickey: 5
topictitle: "Assembly"
---

## Accompanying material

* Slides 2017: [UBC - De novo Assembly 2017](./Assembly2017.pdf)
* Slides 2018: [UBC - De novo Assembly 2018](./Assembly2018.pdf)
* Slides 2019: [UBC - De novo Assembly 2019](./Assembly2019.pdf)
* Background reading: Comparison of the two major classes of assembly algorithms: overlap-layout-consensus and de-bruijn-graph [Paper](./background_reading/Briefings in Functional Genomics-2011-Li-bfgp-elr035.pdf). Briefings in Functional Genetics 2011.
* Sense from sequence reads: methods for alignment and assembly [Paper](./background_reading/Flicek&Birney2009.pdf). Flicek and Birney. Nature methods supplement 2009.
* [Velvet Manual 1.1](./background_reading/Manual.pdf). Daniel Zerbino, 2008.
* Stacks: Building and Genotyping Loci De Novo From Short-Read Sequences [Paper](./background_reading/Stacks.pdf). Catchen, Amores, Hohenlohe, Cresko, Postlethwait. G3 Genes Genomes Genetics, 2011.
* Trinity: reconstructing a full-length transcriptome without a genome from RNA-Seq data [Paper](./background_reading/nihms292662.pdf). Grabherr, Haas, Yassour, Levin, Thompson et al. Nat Biotechnol. 2013.
* The impact of third generation genomic technologies on plant genome assembly [Paper](./Jiao.pdf). Jiao & Schneeberger. Curr. Op. Plant Biol. 2017.


# Topic 5 Assembly

## Code break questions 

1. Write a one liner to find all the overlaps exactly 4 bp in length between CTCTAGGCC and a list of other sequences in the file /home/biol525d/Topic_5/data/overlaps.fa

2. Find all the unique 9mers in a fasta sequence /home/biol525d/Topic_5/data/kmer.fa

This might be a tricky one and there are likely many ways to do this. First try out these commands that might help you complete this task. It might also be helpful to determine how many characters are in the sequence (wc -c).
		
Hints: test out the following commands:

```bash
cut -c1- kmer.fa
```

```bash
cut -c1-4 kmer.fa
```
		
```bash
for num in {1..10}
    do
        echo $num >> file.txt
    done
```
What do these commands do? Can you use commands like this to find all the kmers in the sequence?

3. Sort them and keep the unique kmers

Hint: try sort (look up the options)


## Tutorial 

Your goal for today is to assemble a bacterial genome, as best you can, using the program Velvet

This genome was downloaded from GAGE (Genome Assembly Gold-Standard Evaluations) website (http://gage.cbcb.umd.edu/data/index.html): 

Details for this dataset are as follows: 
⋅⋅*Species: Staphylococcus aureus
⋅⋅*Actual genome size: 2,860,307 bp
⋅⋅*Type: Paired end
⋅⋅*Read number (total - including both reads per pair): 1,294,104
⋅⋅*Read size (each read): 101 bp 
⋅⋅*Insert length (sd): 180 bp (+/-20 bp) 

The data is located in /home/biol525d/Topic_5/data/. One file contains the forward read (frag_1.fastq.gz) and the other file contains the reverse read (frag_2.fastq.gz). Each read has a match from the same fragment in the other file and is in the same order in the matching file.

Step 1. Install Velvet COMPLETED 
(see Manual at http://www.ebi.ac.uk/~zerbino/velvet/Manual.pdf):

Velvet is already installed, but this was the command used to install the program:

sudo apt-get install velvet

(note I also placed the zipped program in the /home/biol525d/Topic_5/scripts folder just in case, but it would have to be installed)

Step 2. Enter this folder and unpack these data: COMPLETED 

gunzip -d frag_*.fastq.gz


### Question 1) Given the above information, what is the expected coverage?

This information will be useful when you are running your genome assemblies.

Step 3. The first program you need to run is velveth

Velveth takes in a number of sequence files, produces a hashtable (i.e. all the kmers), then outputs two files in an output directory (creating it if necessary). These files are called Sequences and Roadmaps, and are necessary to velvetg. 

The syntax to run velveth is as follows:

velveth <output_directory> <hash_length>

To ensure that each k-mer cannot be its own reverse complement, k (i.e. hash length or kmer length) must be odd.

For all the options simply type velveth and refer to the velvet manual for more details.

Here is an example command line:

Make a directory in your home directory for your output

```bash
mkdir ~/Topic_5
```

Move into the ~/Topic_5 directory and run velveth

```bash
velveth sa_assembly21 21 -shortPaired -fastq  -separate /home/biol525d/Topic_5/data/frag_1.fastq /home/biol525d/Topic_5/data/frag_2.fastq
```

21 is the kmer length
-shortParied specifies the types of reads (paired end)
-fastq is the type of sequence format
-separate tells the programs that the paired reads are in two separate files (one with read 1 (frag_1) and one with read 2 (frag_2))

Run the above command.

A Roadmap file is produced. This file is used as input for the next stage of velvet. For each k-mer observed in the set of reads, the hash table records the ID of the first read encountered containing that k-mer and the position of its occurrence within that read. This file rewrites each read as a set of original k-mers combined with overlaps with previously hashed reads. For more information see:

http://homolog.us/blogs/blog/2011/12/06/format-of-velvet-roadmap-file/

This command needs a great deal of memory and may not work on your instance. Therefore, I have run this command using two different kmers (31 and 21). You can use the Roadmaps that I generated complete the rest of the tutorial.

Step 4. The next step is to make make the graph, simplify the graph, correct errors and resolve repeats. This is done by velvetg.

velvetg is run as follows:

velvetg <velveth output_directory>

For all the options simply type velvetg and refer to the velvet manual for more details.

Here is an example command line with no options:

```bash
velvetg sa_assembly21
```

Run the above command.

Now contigs.fa appears in your output directory (sa_assembly21). This file contains your assembled genome. A log file, with information on your input parameters as well as some basic metrics of your assembly, and a stats.txt file with information on each node is present also appear. Note that node lengths are given in k-mers (see below). 

Units in Velvet:

Velvet measures and reports lengths in overlapping k-mers. Although not intuitive at first sight, this unit system allows for consistency throughout Velvet’s output. 

The relationship between coverage and k-mer coverage is defined by the following equation:

Ck = C ∗ (L−k+1)/L

Where C=coverage, L=read length, and k=kmer length. 

Similarly, statistics derived from lengths are also subjected to this transformation. It is as simple as adding k-1 to the reported length to recover the length in bp. If the median coverage of an assembly is reported as Ck read k-mers per contig k-mer it is corresponds in fact to roughly CkL/(L + k −1) read basepair per contig basepair (assuming that contigs are significantly longer than the hash length).

### Question 2) For a k-mer of 21 what is the k-mer coverage for this genome assembly?


Step 5. Assess assembly 

There are several basic metrics to quantify the quality of a genome assembly. N50 is a common statistic similar to a mean or median contig length, but has greater weight given to the longer contigs. It is defined as the contig length at which half the bases in the genome are in contigs that size or larger. Other metrics include the longest contig size, the total size of the assembly and total contig number. 

### Question 3) Can you think of other ways to assess assembly quality? What might be the trouble with only focusing on maximizing N50? Discuss this with your group.

You can also quantify the assembly metrics using the perl script fasta_lengths_workshop.pl (in addition to the output provided from velvet). You will not need to make the adjustment for k-mer length as you would for the velvet assembly statistics.

Run as follows in the sa_assembly21 directory:

```bash
perl  /home/biol525d/Topic_5/scripts/fasta_lengths_workshop.pl contigs.fa
```

This script will produce two files. The first will be called contigs.fa.lengths, which is the length of each contig in the assembly. The second is contigs.fa.stats which provides: 

1) The number of contigs in the assembly file contigs.fa
2) The average length of the contigs
3) The total number of bp in the assembly file
4) The median contig length
5) The minimum contig length
6) The maximum contig length
7) The N50

### Question 4) Quantify the assembly metrics for your first assembly that you ran without any options. In your group of four, each person should pick different sets of parameters to run. Compare the resulting assemblies with one another and discuss which ones seemed to have improved the assembly and why that might be. Be prepared to share your findings with the class. 

Be careful not to use all your hard disk space or RAM. Velvet is relatively fast, but a memory hog. When you are finished the tutorial, you may want to remove files produced by Velvet from your drive to make space. 


### Recap: 

You can check how much space is on your drive by typing: 

```bash
df -h
```

This will provide you information on how large the directory you are in is:

```bash
du -sh 
```

This will tell you how much memory is being used currently by the server:

```bash
free -g 
```

To remove files you can type: 

```bash
rm <filename>
```

To remove all files in your current directory enter the directory and type the following (Be careful not to delete something you need! Every file in that directory will be GONE!)

```bash
rm *
```

To remove an empty directory type:

```bash
rmdir <directory name>
```

To move files or directories:

```bash
mv <old file name> <new file name>
```

### Back to the tutorial

If you want to modify the k-mer you must run velveth again and replace 21 with a new number. Velvet only allows k-mers up to 31 bp in length. 

Typically, the longer the k-mer, the better the assembly, until you hit the point of too little coverage. Discuss why this might be with your group.

If you want to modify other parameters you can simply run velvetg without rerunning velveth. (Note that if you re-run fasta_lengths_workshop.pl and velvetg using the same filenames in the same directory the files will be written over, so either record your run options and assembly metrics, rename the output files or move them to a new directory to retain the information).

Some potential parameters to modify include (see the manual for details):
-min_contig_lgth
-cov_cutoff
-ins_length 
-ins_length_sd 
-exp_cov (note that for this parameter you can include an estimated expected k-mer coverage or ask velvet to estimate it from the data by typing -exp_cov auto)

Be sure to test -exp_cov. For a more detailed explanation of this parameter see:

http://homolog.us/blogs/blog/2012/06/08/an-explanation-of-velvet-parameter-exp_cov/
