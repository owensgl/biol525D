#install the package edgeR
source ("http://www.bioconductor.org/biocLite.R")

biocLite("edgeR")

#find your working directory:
getwd()

#paste it in here (i.e. replace my path with yours):
setwd ("/Users/kayhodgins/Dropbox/Documents/bioinformatics_workshop/bioinformatics_examples/data_for_exercises/rnaseq")

#load the libraries you will need 
library ("edgeR")
library ("gplots")

#read in the data
data <- read.table ("cold_hot_expression.txt")

#extract and turn the column names into a factor
treat <- as.factor (sapply (strsplit(colnames(data),split = ""),"[[",1))

#make a DGElist object out of the data
list1 <- DGEList (counts = data, group = treat)

#calculate the normalization factors to adjust the effective library size relative to other libraries in the dataset
list1 <- calcNormFactors (list1)

#calculate the counts per million
cpm.list1 <- cpm(list1)

#filter out the genes with less than 1 cpm in 6 or fewer libraries (a somewhat arbitrary choice). Genes are usually dropped if they can't possibly be expressed in all the samples for any of the conditions.
list2 <- list1[rowSums(cpm.list1 > 1) >= 6,]
cpm.list2 <- cpm.list1[rowSums(cpm.list1 > 1) >= 6,]

#generate a multi-dimensional scaling plot
col_treat <- as.character (treat)
col_treat [col_treat == "C"] <- "blue"
col_treat [col_treat == "H"] <- "red"

#plot the MDS graph
plotMDS (list2, col = col_treat)

#Set the model to use. This one includes the intercept, but other models can be specified that omit the intercept or that have more complex designs. See EdgeR manual for details.
design <- model.matrix (~treat)

#fit the common + tagwise dispersion models
list2 <- estimateGLMCommonDisp (list2, design)
list2 <- estimateGLMTagwiseDisp (list2, design)

#fit a GLM to the data using the tagwise dispersion estimate
glm.list2 <- glmFit (list2, design, dispersion = list2$tagwise.dispersion)
lrt.list2 <- glmLRT (glm.list2)

#get the topTags out of the model
top <- topTags (lrt.list2, n = 1000)$table

fdr<-p.adjust(lrt.list2$table$PValue, method='fdr')
hist(fdr,  breaks=1000)
dim(lrt.list2$table[fdr<0.05,])

#make a heatmap by getting the counts per million from each gene and turning them relative proportions (columns add up to 1)
sub1 <- colSums (cpm.list2)
sub2 <- matrix (rep(sub1,nrow (cpm.list2)), c (nrow (cpm.list2),ncol(cpm.list2)),byrow = TRUE)
sub3 <- cpm.list2 / sub2

#subset this matrix to just get the genes from above
names1 <- row.names (sub3)
names2 <- row.names (top)
index2 <- names1 %in% names2
heatmap1 <- sub3[index2,]

#play around with the options to make the plot fit what you like for options type ?heatmap.2
heatmap.2 (heatmap1, trace = "none", scale = "row", Colv = FALSE, labCol = treat,cexCol = 1, dendrogram = "none", labRow = FALSE, colsep =c(6), sepcolor = "white", sepwidth = 0.03)


#plot the individual expression values from a single gene:
subdat1 <- data["comp520_c0",]
stripchart (log10(as.numeric (subdat1)) ~ treat, vertical = T, xlim = c (0.5,2.5), method = "jitter", ylab = "Log10 Expression count")


