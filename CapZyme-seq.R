# CapZyme-seq data analysis pipeline R part (Sherwood et al, 2021)
# Jeppe Vinther, 2021-07-29

#### Requirements
# R (https://www.r-project.org/)
# featureCounts output file from a CapZyme experiment
# DESeq2 package (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# tidyverse (https://www.tidyverse.org/)
####


# Preparation for CapZyme-seq analysis ---- 
library(DESeq2)
library(tidyverse)
library(reshape2)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Example data
# Example data is 10 fastq files (reduced in size) from CapZyme-seq libraries prepared from
# J6/JFH1 infected Huh-7.5 cells and treated with either AtNUDX23 or Control (no treatment) 
# AtNUDX23: 02, 05, 08, 11, 14
# Control: 03, 06, 09, 12, 15
# CapZyme-seq.sh file shows how to produce the featureCounts file from fastq files
# Example count file is located $path/FeatureCount_counts.txt

# Read count file, precomputed available at erda.ku.dk
counts <- read_delim("https://erda.ku.dk/archives/e63276ecfaa314afaf8305d8e09fbe95/FeatureCount_counts.txt", 
                     "\t", escape_double = FALSE, comment = "#", 
                     trim_ws = TRUE)


# DESeq2 paired analysis----
# Prepare count matrix
countData <- as.matrix(counts[,c("02.sam", "05.sam", "08.sam", "11.sam", "14.sam",
                                 "03.sam", "06.sam", "09.sam", "12.sam", "15.sam")])
rownames(countData) <- counts$Geneid
head(countData)

# Prepare column-data dataframe (information about samples)
condition <- c("AtNUDX23","AtNUDX23","AtNUDX23","AtNUDX23","AtNUDX23",
               "Control","Control","Control","Control","Control")
sample <- c("Sample1","Sample2","Sample3","Sample4","Sample5","Sample1","Sample2","Sample3","Sample4","Sample5")

colData <- data.frame(condition,sample)
rownames(colData) <- c("02.sam", "05.sam", "08.sam", "11.sam", "14.sam",
                       "03.sam", "06.sam", "09.sam", "12.sam", "15.sam")
colData$condition <- factor(colData$condition)
colData$sample <- factor(colData$sample)
colData$condition <- relevel(colData$condition,"Control")

# The columns of the count matrix and the rows of the column-data dataframe should be in the same order
all(rownames(colData) == colnames(countData))

# Make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ sample + condition)


# Run DESeq analysis
dds <- DESeq(dds)


# Make result for the AtNUDX23/Control contrast
res <- results(dds, contrast=c("condition","AtNUDX23","Control"))
summary(res)



resultsAtNUDX23 <- as.data.frame(res[order(res$pvalue),])
resultsAtNUDX23$Geneid <- row.names(resultsAtNUDX23)

# Write result to table----
resultsAtNUDX23 <- as.data.frame(res[order(res$pvalue),])
resultsAtNUDX23$Geneid <- row.names(resultsAtNUDX23)
setwd("INSERT_PATH_TO_FOLDER")
setwd("/binf-isilon/vintherlab/jvinther/scratch")
write.table(resultsAtNUDX23, file = "CapZyme_AtNUDX23_results.txt", sep = "\t", row.names = F, dec = "," )

# Plot results----

# Standard DESeq2 plots
plotDispEsts(dds)
plotMA(res, ylim=c(-2.5,2.5), main = "AtNUDX2 enrichment")
plotCounts(dds = dds, gene = "HCV_J6-JFH1-c2-plus", intgroup="condition")
plotCounts(dds = dds, gene = "HCV_J6-JFH1-c2-minus", intgroup="condition")


# Make MD plot
plot(log(resultsAtNUDX23$baseMean),resultsAtNUDX23$log2FoldChange, xlab ="log(mean expression)", 
     ylab= "log2fold(AtNUDX23 counts/Control counts)", main ="AtNUDX23 enrichment",
     ylim=c(-4.5,4.5),col = alpha(colour = cbPalette[1], alpha = .3), pch = 20)
points(log(resultsAtNUDX23$baseMean[grepl(pattern = "ERCC", x = resultsAtNUDX23$Geneid)]),
       resultsAtNUDX23$log2FoldChange[grepl(pattern = "ERCC", x = resultsAtNUDX23$Geneid)], col = "black", bg =cbPalette[3], pch = 21)
points(log(resultsAtNUDX23$baseMean[grepl(pattern = "HCV", x = resultsAtNUDX23$Geneid)]),
       resultsAtNUDX23$log2FoldChange[grepl(pattern = "HCV", x = resultsAtNUDX23$Geneid)], col = "black", bg =cbPalette[5], pch = 21)
points(log(resultsAtNUDX23$baseMean[grepl(pattern = "5S", x = resultsAtNUDX23$Geneid)]),
       resultsAtNUDX23$log2FoldChange[grepl(pattern = "5S", x = resultsAtNUDX23$Geneid)], col = "black", bg =cbPalette[6], pch = 21)
points(log(resultsAtNUDX23$baseMean[grepl(pattern = "BCYRN1", x = resultsAtNUDX23$Geneid)]),
       resultsAtNUDX23$log2FoldChange[grepl(pattern = "BCYRN1", x = resultsAtNUDX23$Geneid)], col = "black", bg =cbPalette[6], pch = 21)
text(log(resultsAtNUDX23$baseMean[grepl(pattern = "HCV_J6-JFH1-c2-plus", x = resultsAtNUDX23$Geneid)]),
     resultsAtNUDX23$log2FoldChange[grepl(pattern = "HCV_J6-JFH1-c2-plus", x = resultsAtNUDX23$Geneid)], labels="HCV+", cex=0.7, font=1, pos=1)
text(log(resultsAtNUDX23$baseMean[grepl(pattern = "HCV_J6-JFH1-c2-minus", x = resultsAtNUDX23$Geneid)]),
     resultsAtNUDX23$log2FoldChange[grepl(pattern = "HCV_J6-JFH1-c2-minus", x = resultsAtNUDX23$Geneid)], labels="HCV-", cex=0.7, font=1, pos=1)
legend("topright",cex = .7, pt.cex= 1, legend = c("HCV","ERCC","polIII"), col = rep("black",3), pt.bg = c(cbPalette[5],cbPalette[3],cbPalette[6]), pch = c(21,21,21))


# Make vulcano plot
plot(resultsAtNUDX23$log2FoldChange, -log10(resultsAtNUDX23$pvalue), ylab ="-log10(p-value)", 
     xlab= "log2fold(AtNUDX23 counts/control counts)", main ="AtNUDX23 enrichment",
     col = alpha(colour = cbPalette[1], alpha = .3), pch = 20)
points(resultsAtNUDX23$log2FoldChange[grepl(pattern = "HCV", x = resultsAtNUDX23$Geneid)], 
       -log10(resultsAtNUDX23$pvalue[grepl(pattern = "HCV", x = resultsAtNUDX23$Geneid)]), col = "black", bg =cbPalette[5], pch = 21)
points(resultsAtNUDX23$log2FoldChange[grepl(pattern = "5S", x = resultsAtNUDX23$Geneid)], 
       -log10(resultsAtNUDX23$pvalue[grepl(pattern = "5S", x = resultsAtNUDX23$Geneid)]), col = "black", bg =cbPalette[6], pch = 21)
legend("topleft",cex = .7, pt.cex= 1, legend = c("HCV","5S rRNA"), col = rep("black",2), pt.bg = c(cbPalette[5],cbPalette[6]), pch = c(21,21))


# Accessory analysis: Analysis of 5 termini sequence ----

#### Requirements
# R (https://www.r-project.org/)
# Wordcount output file from a CapZyme experiment (as described in CapZyme.sh)
# Example data is in $path/data/termini_nt/output
# Reshape2 package (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# stringr package
# tidyverse (https://www.tidyverse.org/)
####
library(tidyverse)
library(reshape2)
library(stringr)

# Preparation for termini nt-seq analysis ---- 
# set path
outputPath <- SET_PATH_TO_OUTPUT_FILE_DIRECTORY
setwd(outputPath)
files <- list.files(outputPath)

# Read termini nt counts
terminiCounts <- data.frame(nt = c("A","G","C","U"))
for(i in files) {
  terminiCounts[i] <- as.numeric(scan(i, what = "numeric"))
}

# Prepare for plotting with ggplot2 and annotate samples
terminiCounts <- melt(terminiCounts, id.vars=c("nt"))
terminiCounts$sample <- as.factor(str_split(terminiCounts$variable, pattern = "_" , n = Inf, simplify = TRUE)[,1])
terminiCounts$treatment <- NA
terminiCounts$treatment[terminiCounts$sample  %in% c("02","05","08","11","14")] <- "AtNUDX23"
terminiCounts$treatment[terminiCounts$sample  %in% c("03","06","09","12","15")] <- "Control"
terminiCounts$treatment <- as.factor(terminiCounts$treatment)


# Annotate replicates
levels(terminiCounts$sample)[levels(terminiCounts$sample) == "02" | levels(terminiCounts$sample) == "03"] <- "rep1"
levels(terminiCounts$sample)[levels(terminiCounts$sample) == "05" | levels(terminiCounts$sample) == "06"] <- "rep2"
levels(terminiCounts$sample)[levels(terminiCounts$sample) == "08" | levels(terminiCounts$sample) == "09"] <- "rep3"
levels(terminiCounts$sample)[levels(terminiCounts$sample) == "11" | levels(terminiCounts$sample) == "12"] <- "rep4"
levels(terminiCounts$sample)[levels(terminiCounts$sample) == "14" | levels(terminiCounts$sample) == "15"] <- "rep5"
names(terminiCounts)[names(terminiCounts)=="value"] <- "Percentage"

# Plot nt termini counts
ggplot(terminiCounts, aes(fill=nt, y=Percentage, x=sample)) + 
  geom_bar(position="fill", stat="identity")  +
  scale_fill_manual(values = cbPalette[c(5,3,8,7)]) +
  ggtitle("Start nucleotide J6/JFH1+")
  
# Accessory analysis: Analysis of sequencing depth ----

#### Requirements
# R (https://www.r-project.org/)
# Depth output file from a CapZyme experiment (example data: $path/data/featureCounts/depth_file.txt)
# featureCounts summary file from a CapZyme experiment (example data: $path/featureCount_counts.txt.summary)
# DESeq2 package (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# tidyverse (https://www.tidyverse.org/)
####
# Function to normalize the sequencing depth to per 10^6 assigned reads
normaliseDepth <- function(count_sum_file, depthFile, ...) {
  depthNorm <- data.frame(chr = depthFile$`#CHROM`, POS = depthFile$POS)
  depthNorm[, colnames(depthFile)[3:ncol(depthFile)]] <- depthFile[,colnames(depthFile)[3:ncol(depthFile)]]*10E6/as.data.frame(lapply(counts_sum[1,gsub(pattern = ".bam", replacement = ".sam",x = colnames(depthFile)[3:ncol(depthFile)])],rep,nrow(depthFile)))
  depthNorm
}
# Function to plot sequencing depth for specific RNA
plotDepth <- function(RNA_ID, BAM_FILE, depthFile, ...) {
  depthFile <- depthFile[depthFile$chr == RNA_ID,]
  o <- order(depthFile$POS[depthFile$chr == RNA_ID] )
  depthFile <- depthFile[o,]
  plot(as.vector(depthFile$POS),depthFile[, (names(depthFile) %in% BAM_FILE)],
       type="h",
       main = paste(RNA_ID, BAM_FILE, sep = "\n"),
       xlab = "Position", 
       ylab = "Sequencing depth",
       ...)
  
}
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


countsSum <- read.delim("/binf-isilon/vintherlab/jvinther/scratch/featureCount_counts.txt.summary", check.names = F,header = TRUE,sep = "\t")
depthFile <- read.delim("/binf-isilon/vintherlab/jvinther/scratch/data/featureCounts/depth_file.txt", check.names = F,header = TRUE,sep = "\t" )

# Normalise the sequencing depth to library depth (10^6 assigned reads)
depthNorm <- normaliseDepth(countsSum,depthFile)


# Plotting sequencing depth
par(mfrow=c(2,2))

plotDepth(RNA_ID = "HCV_J6-JFH1-c2-plus",
           BAM_FILE = "05.bam", 
           depthFile = depthNorm, 
           col = cbPalette[6], 
           xlim= c(1,9678), 
           ylim=c(1,3500),
           bty="n",
           lwd=1) 


plotDepth(RNA_ID = "HCV_J6-JFH1-c2-minus",
           BAM_FILE = "05.bam", 
           depthFile = depthNorm, 
           col = cbPalette[4], 
           xlim= rev(c(1,9678)), 
           ylim=rev(c(1,600)),
           bty="n",
           lwd=1,
           axes = TRUE) 

plotDepth(RNA_ID = "HCV_J6-JFH1-c2-plus",
           BAM_FILE = "06.bam", 
           depthFile = depthNorm, 
           col = cbPalette[6], 
           xlim= c(1,9678), 
           ylim=c(1,3500),
           bty="n",
           lwd=1) 

plotDepth(RNA_ID = "HCV_J6-JFH1-c2-minus",
          BAM_FILE = "06.bam", 
          depthFile = depthNorm, 
          col = cbPalette[4], 
          xlim= rev(c(1,9678)), 
          ylim=rev(c(1,600)),
          bty="n",
          lwd=1,
          axes = TRUE) 

