# CapZyme data analysis
This page contains information related to the analysis of the data in _Sherwood et al, 2021_. In this paper, we use the CapZyme-seq method to identify RNAs that are 5' capped with the metabolic cofactor Flavine Adenine Dinucleotide (FAD).

The CapZyme experiment compares a Control sample to a sample treated with an enzyme (AtNUDX23, Rpp or other), which specifically generate mono-phosphate at the 5’ terminal of RNAs depending on their 5' termini. Subsequent 5’ monophosphate specific ligation of a sequencing adaptor to RNAs allows the RNAs with a 5' modification (such as FAD) to be enriched in a sequencing experiment. This page decsribes the method used in _Sherwood et al, 2021_ to analyse the data.



![CapZyme-seq](https://user-images.githubusercontent.com/42373129/126970748-336ac8fd-cec4-4ba5-b63a-72f0b5e56ac3.png)

<img width="400" alt="MD plot" src="https://user-images.githubusercontent.com/42373129/126970748-336ac8fd-cec4-4ba5-b63a-72f0b5e56ac3.png">

CapZyme-seq bash data analysis, described in CapZyme-seq.sh (https://github.com/jeppevinther/CapZyme/blob/main/CapZyme-seq.sh):

1. Trimming with cutadapt
2. Pseudo mapping to entire transcriptome with kallisto
3. Mapping to expressed mRNAs + viral RNAs + small RNAs with Bowtie2
4. BAM to SAM files with samtools
5. Counting 5' termini reads with FeatureCount

CapZyme-seq R data analysis, described in CapZyme-seq.R (https://github.com/jeppevinther/CapZyme/blob/main/CapZyme-seq.R)

6. Importing count files into R
7. Analysis with DeSeq2
8. Plotting the data

Using the scaled-down example data, the main finding of _Sherwood et al, 2021_ can be reproduced:
<img width="400" alt="MD plot" src="https://user-images.githubusercontent.com/42373129/127151633-de031479-0f8e-4624-b9f8-92cbc7f3165f.png">

Other analyses, described in XXXXXXX.sh:

1. Making sequencing depth files
2. Analysis of 5' termini sequence
3. Making wig files for uploading to UCSC as custom tracks

