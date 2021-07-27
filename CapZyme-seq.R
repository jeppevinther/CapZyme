# CapZyme-seq data analysis pipeline (Sherwood et al, 2021)

#### Requirements
# R (https://www.r-project.org/)
# featureCounts output file from a CapZyme experiment
# DESeq2 package (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
# tidyverse (https://www.tidyverse.org/)
####


# Preparation for analysis ---- 
library(DESeq2)
library(tidyverse)
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Example data
# Example data is 10 fastq files (reduced in size) from CapZyme-seq libraries prepared from
# J6/JFH1 infected Huh-7.5 cells and treated with either AtNUDX23 or Control (no treatment) 
# AtNUDX23: 02, 05, 08, 11, 14
# Control: 03, 06, 09, 12, 15
# CapZyme-seq.sh file shows how to produce the featureCounts file from fastq files

# Read count file, 
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

