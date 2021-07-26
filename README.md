# CapZyme data analysis
This page contains information related to the analysis of the data in _Sherwood et al, 2021_. In this paper, we use the CapZyme-seq method to identify RNAs that are 5' capped with the metabolic cofactor Flavine Adenine Dinucleotide (FAD).

The CapZyme experiment compares a Control sample to a sample treated with an enzyme (AtNUDX23, Rpp or other), which specifically generate mono-phosphate at the 5’ terminal of RNAs depending on their 5' termini. Subsequent 5’ monophosphate specific ligation of a sequencing adaptor to RNAs allows the RNAs with a 5' modification (such as FAD) to be enriched in a sequencing experiment. This page decsribed the method used in _Sherwood et al, 2021_ to analyse the data. 

![CapZyme-seq](https://user-images.githubusercontent.com/42373129/126970748-336ac8fd-cec4-4ba5-b63a-72f0b5e56ac3.png)


The flow of the CapZyme-seq data analysis described in XXXXXXX.sh:
1. Trimming with cutadapt
2. Pseudo mapping to entire transcriptome with kallisto
3. Mapping to expressed mRNAs + viral RNAs + small RNAs with Bowtie2
4. Counting 5' termini reads with FeatureCount
5. Importing count files into R
6. Analysis with DeSeq2.

