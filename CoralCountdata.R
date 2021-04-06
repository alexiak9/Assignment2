#author: Alexia Kotorov
#step 1: set the working directory
setwd("~/Users/lexikotorov/Desktop/School yeehaw/Assignment 2/") 

#step 2: Install all necessary packages for code. 
#These packages contain all the necessary tools for DESeq2 analysis of the imported data files. These packages are standard for such workflow.
library(DESeq2) 
library(affycoretools)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library("dplyr")

#step 3: import the data and read the counts in it 
#The data is from the ex situ experiments. Ex situ experiments are an important part of ecological studies as they reduce the variation present in in-situ experiments.
#The data in the countsfile.flow file is of the 4 high flow tanks, and 3 low flow tanks. The file contains transcriptomic data.
countData <- read.table("countsfile.flow")
head(countData)
length(countData[,1])
#there are 27807 counts
#from the header can see that out of the four categories, the data table includes high.flow and low.flow

#step 4: rename the row names to restructure the data table, so that it is easier to follow.
row.names(countData)=sub("TRINITY_GG_", "isogroup", rownames(countData))
head(countData)

#step 5: create a folder in which to carry out outlier functions
#set a variable (V) to the folder path so can access in functions
#It is crucial to remove any outliers prior to the rest of the analysis and to keep records of the outlier data in a separate file. 
#Outliers can be caused due to RNA contamination or technical errors during RNA processing in building the transcriptomic profiles, and would impact our subsequent analysis. 
dir.create("/Users/lexikotorov/Desktop/School yeehaw/Assignment2/outlier")
setwd("/Users/lexikotorov/Desktop/School yeehaw/Assignment2/outlier")
v=setwd("/Users/lexikotorov/Desktop/School yeehaw/Assignment2/outlier")

newcountData=round(countData, 0)
treat=c("High.Flow", "High.Flow", "High.Flow", "High.Flow", "Low.Flow", "Low.Flow", "Low.Flow")
g=data.frame(treat)
g
colData= g

dds=DESeqDataSetFromMatrix(countData=newcountData,
                           colData = g,
                           design = ~treat)

vsd.ge=assay(vst(dds))
rl=vst(dds)
e=ExpressionSet(assay(rl), AnnotatedDataFrame(as.data.frame(colData(rl))))
arrayQualityMetrics(e,outdir=v,intgroup=c("treat"),force=T)
dev.off()

#restart R after removing outliers once 
#outlier report has been saved in outlier folder 

#step 6: re-read countdata to see what outliers were removed and switch back to wd
#This step is important again to keep track of what data we are actually analyzing. In working with such large data sets/data files, it is important to remove outliers, know how many outliers were removed, and how many data points were maintained. 
#After filtering the data points for outliers, we are now ready to officially begin analyzing.
setwd("~/Users/lexikotorov/Desktop/School yeehaw/Assignment 2/")
library("DESeq2")
library("ggplot2")
library("dplyr")
BiocManager::install("tximport")
library(tximport)

countData <- read.table("countsfile.flow")
head(countData)
length(countData[,1])

row.names(countData)=sub("TRINITY_GG_", "isogroup", rownames(countData))
head(countData)

#Step 7: barplot the values you have, the total counts per each column or sample group
#This barplot of the values lets us see how many counts we have per each sample group. The high flow derived reads are indicated in coral on the barplot, and the low flow reads are indicated in blue.
#We can see from the bar plot that the most counts came from the High.Flow_1 group. Variation in count reads ___
totalCounts=colSums(countData)
totalCounts
barplot(totalCounts, col=c("coral","coral","coral","coral","blue","blue","blue"), ylab="raw counts")

min(totalCounts) #22023.1
max(totalCounts) #3709559

#step 8: DESeq analysis of differential expression
#DESeq analysis shows how genes are differentially expressed based on our variable of choice (in this case: high vs. low flow).
#The results summary shows the log fold change and shows results for each row, or each sample. 
#I had to assign a new newcountData variable that rounded up decimal points within the data so that I was only inputing integers. While this is not the most accurate way to perform this analysis, rounding up via this code & using only integers was sufficient for the sake of this experiment. 
newcountData=round(countData, 0)
treat=c("High.Flow", "High.Flow", "High.Flow", "High.Flow", "Low.Flow", "Low.Flow", "Low.Flow")
g=data.frame(treat)
g
colData<- g

dds<-DESeqDataSetFromMatrix(countData=newcountData, colData=colData, design=~treat)

#one step DESeq
dds<-DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing

head(dds)
res<- results(dds)

#step 9: dispersion plot to look at variance of counts 
#Dispersion plots provide a visual for how much read counts of a specific RNA sample differ from the mean.
#This visual is important in starting to see how the treatment of interest's genetic profile (high flow in this case) differs from the baseline (low flow)
plotDispEsts(dds, main="Dispersion plot of Coral in Different Flow Microenvironments")

#step 10: high flow vs. low flow pairwise comparison
colData$HighFlow<-factor(colData$treat, levels=c("High.Flow","Low.Flow"))
##second term is the "control", in this case being low flow as we want to see how high flow differs from low flow
resHighFlow <- results(dds, contrast=c("treat","High.Flow","Low.Flow"))

#Step 11: We want to see how many FDR < 10%?
#The false discovery rate is set to a certain threshold to ensure that the genes deemed differentially upregulated or downregulated are actually due to changes in gene expression profile.
#Without adjusting for the FDR, we might include "false positives" or genes that aren't actually upregulated or of interest to us. 
table(resHighFlow$padj<0.01)
summary(resHighFlow)
nrow(resHighFlow[resHighFlow$padj<0.05 & !is.na(resHighFlow$padj),])  
# The above code tells us num significantly differentially expressed genes excluding the no/low count genes

#Step 12: plotMA 
#This function gives us the scatter plot of the log fold change between the data points. With this function, we get a visual of the comparison of expression profiles, or the comparison of read outs.
#Taking the adjusted pairwise comparison between the high and low flow rate data (resHighFlow) and plugging that into the plotMA function gives us this data visual.
#We can also read out, as is shown in line 135 and 136, how many reads were upregulated or downregulated, therefore extending beyond just a visual plot & providing us with concrete numerical data.
#The output data table is stored in the variable results and the summary file is output as "High Flow.txt"
dev.off()
plotMA(resHighFlow, main="Low Flow vs. High Flow")
plotMA(resHighFlow, main="Low Flow vs. High Flow", ylim=c(-2,2))

results <- as.data.frame(resHighFlow)
head(results)

nrow(resHighFlow[resHighFlow$padj<0.1 & resHighFlow$log2FoldChange > 0 & !is.na(resHighFlow$padj),])
nrow(resHighFlow[resHighFlow$padj<0.1 & resHighFlow$log2FoldChange < 0 & !is.na(resHighFlow$padj),])
#UP 9237
#DOWN 360

write.table(resHighFlow, file="High Flow.txt", quote=F, sep="\t")

cd <- read.table("High Flow.txt")
head(cd)

#step 11: make GO table for the genes differentially expressed in high flow samples, as determined from the dispersion plot and plotMA
#This function gives us and saves within the go_input_High.Flow table the GO enrichment rankings. 
head(cd)

library(dplyr)
cd
go_input_High.Flow = cd %>%
  tibble::rownames_to_column(var = "iso") %>%
  mutate(mutated_p = -log(pvalue)) %>%
  mutate(mutated_p_updown = ifelse(log2FoldChange < 0, mutated_p*-1, mutated_p*1)) %>%
  na.omit() %>%
  select(iso, mutated_p_updown)

#Step 12: p value analysis for the analysis we completed above for high flow
valHF=cbind(resHighFlow$pvalue, resHighFlow$padj)
head(valHF)
colnames(valHF)=c("pval.HighFlow", "padj.HighFlow")
length(valHF[,1])
table(complete.cases(valHF))

#Step 13: making an rlogdata and pvals table 
rlog=rlogTransformation(dds, blind=TRUE) 
rld=assay(rlog)
head(rld)
colnames(rld)=paste(colData$treat)
head(rld)
length(rld[,1])

rldpvals=cbind(rld,valHF)
head(rldpvals)
dim(rldpvals)
table(complete.cases(rldpvals))


write.csv(rldpvals, "HighFlow_RLDandPVALS.csv", quote=F)

colnames(rld)=paste(colData$treat)
head(rld)

#step 14: sample distance heatmap
#Heatmaps are another incredibly useful tool for visualizing large data sets and comparing treatment types. 
#This heatmap compares again the two treatment types, generating a matrix for visualizing similarity or dissimilarity between profiles of sample types.
library(RColorBrewer)
sampleDists <- as.matrix(dist(t(rld)))
library(gplots)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          margin=c(10, 10), main="Comparing FLow Type Treatment on Transcriptome Profiles")


