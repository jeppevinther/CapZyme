# Bash CapZyme-seq data analysis pipeline (Sherwood et al, 2021)

#### Requirements
# cutadapt (added to the path) https://cutadapt.readthedocs.io/en/stable/index.html
# kallisto (added to the path) https://pachterlab.github.io/kallisto/download.html
# samtools (added to the path) http://www.htslib.org/download/
# bowtie2 (added to the path) http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# awk
# Fastq files from CapZyme experiment
# Fasta file with the sequences of the RNAs of interest
####


##############
# Preparation for analysis 
# Make directory for analysis
mkdir directory
cd /directory

#Put the path to the directory in $path
path=PATH_TO_THE_DIRECTORY

cd $path
mkdir data
cd data
# Copy your fastq files to the data directory

# Assign the sample indexes to the $samples variable
# Example data is 10 fastq files (reduced in size) from CapZyme-seq libraries prepared from
# J6/JFH1 infected Huh-7.5 cells and treated with either AtNUDX23 or Control (no treatment) 
# AtNUDX23: 02, 05, 08, 11, 14
# Control: 03, 06, 09, 12, 15
samples='02 03 05 06 08 09 11 12 14 15'

# To download example fastq files (http://doi.org/10.17894/ucph.d4789e55-479c-4f82-8095-49b1af68ed5a)
for lab_number in $samples
do
wget -nH --no-parent -e robots=off https://erda.ku.dk/archives/9a23ff7d54511e58ebf940702f2c2cd0/"$lab_number".fastq.gz
done
wait


# Optional
# If necessary concatenate fastq files from the same index into one fastq file
for i in 1 2 3 4
do
zcat "$i"_file1.fastq.gz "$i"_file2.fastq.gz > "$i".fastq.gz
done
wait


##############
# Trimming with cutadapt 
# Remove adapters and reads shorter than 15 (depends on the method used for library preparation,
# here standard Illumina adapter)
# Filter reads for quality, for NextSeq sequencing use --nextseq-trim=20

for lab_number in $samples
do
cd $path/data
mkdir $lab_number
cd $lab_number
nice cutadapt -a AGATCGGAAGAGCACACGTCT -q 20 -m 15 $directory/A0"$lab_number"*.fastq.gz -o "$lab_number"_trimmed.fastq.gz 1> "$lab_number"_cutadapt.error &
done
wait

# Collect number of reads and trimming percentage in cutadapt_trimming_info.txt
cd $path
touch cutadapt_trimming_info.txt
for lab_number in $samples
do
cd $path/data/$lab_number
echo $lab_number >>$path/cutadapt_trimming_info.txt
grep -F  'Total reads processed:' *_cutadapt.error >>$path/cutadapt_trimming_info.txt
grep -F  'Reads with adapters:' *_cutadapt.error >>$path/cutadapt_trimming_info.txt
grep -F  'Reads written (passing filters):' *_cutadapt.error >>$path/cutadapt_trimming_info.txt
done
wait



##############
# Pseudomapping with kallisto 
# This step is to avoid mapping to sequences of mRNAs that are not expressed
# The fastq are combined and aligned to a kallisto index

# Concatenate fastqs
cd $path/data
cat */*trimmed.fastq.gz > all_fastq_files.fastq.gz

# Download kallisto files (curtesy of the Pachter lab)
cd $path/data
mkdir sequences
cd sequences
wget https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/ensembl-96/homo_sapiens.tar.gz
tar -xvf homo_sapiens.tar.gz

# Pseudomapping
kallisto_index=$path/data/sequences/homo_sapiens/transcriptome.idx
kallisto_fasta=$path/data/sequences/homo_sapiens/Homo_sapiens.GRCh38.cdna.all.fa

cd $path/data
kallisto quant -i $kallisto_index --single -l 200.0 -s 30.0 --threads=6 -o kallisto all_fastq_files.fastq.gz

# Extract mRNAs expressed TPM > 2 and make fasta file for these RNAs
cd $path/data/kallisto
awk '$5>2' abundance.tsv > expressed_transcripts.txt #filter for RNAs with TPM > 2
awk -F'\t' '  NR > 1 { print $1 }' expressed_transcripts.txt > IDs.txt #get IDs
cut -c 1- IDs.txt | xargs -n 1 samtools faidx $kallisto_fasta > expressed.fa #write fasta


##############
# Prepare fastafile for mapping
# Combine your fasta sequences of interest in on fasta file 

# optional download sample fasta files (http://doi.org/10.17894/ucph.d4789e55-479c-4f82-8095-49b1af68ed5a)
cd $path/data/sequences
for file in ERCC92.fa \
            hg19-tRNAs_nonredundant_nointron_CCA_hisG.fa \
            human_rrna.fasta \
            human_srna_GRCh38_p10_uniqed.fa \
            mir_hairpin_human.fa \
            J6-JFH1-clone2.fasta
do
wget -nH --no-parent -e robots=off https://erda.ku.dk/archives/9a23ff7d54511e58ebf940702f2c2cd0/"$file"
done
wait

# Make combined fasta file for mapping
cd $path/data/sequences
cat J6-JFH1-clone2.fasta \
    human_rrna.fasta \
    human_srna_GRCh38_p10_uniqed.fa \
    mir_hairpin_human.fa \
    hg19-tRNAs_nonredundant_nointron_CCA_hisG.fa \
    ERCC92.fa \
    $path/data/kallisto/expressed.fa \
    > RNAs.fa


##############
# Mapping reads with bowtie2


# Prepare bowtie index for mapping
cd $path/data/sequences
nice bowtie2-build RNAs.fa bowtie2_index

# Mapping to the forward strand using --very-sensitive setting and end-to-end.
for lab_number in $samples
do
cd $path/data/"$lab_number"
bowtie2 -p32 --norc --very-sensitive -x $path/data/sequences/bowtie2_index -U "$lab_number"_trimmed.fastq.gz 2>bowtie2.error | gzip > "$lab_number"_mapped.sam.gz
done
wait

# Collect mapping statistics in bowtie2_mapping_info.txt
cd $path
touch bowtie2_mapping_info.txt
for lab_number in $samples
do
cd $path/data/$lab_number
echo $lab_number >> $path/bowtie2_mapping_info.txt
cat bowtie2*.error >> $path/bowtie2_mapping_info.txt
done
wait


##############
# Counting reads with featureCounts

# A SAF file is prepared to use for input into FeatureCount (defining regions for counting)
# Make SAF file based on the the combined fasta file
cd $path/data/sequences
samtools faidx RNAs.fa
awk -F'\t' '{ print $1,$1,1,$2,"+" }' OFS='\t' RNAs.fa.fai > RNAs.saf #Full length of RNAs
awk -F'\t' '$4>100{$4="100"}1' OFS='\t' RNAs.saf > RNAs_100.saf # First 100 nts of RNAs


# Make sam files for featurecount analysis
cd $path/data
mkdir featureCounts
for lab_number in $samples
do
cd $path/data/"$lab_number"
nice gunzip -c "$lab_number"_mapped.sam.gz > $path/data/featureCounts/"$lab_number".sam
done
wait

# Count reads mapping within the first 100 nt of each RNA using featureCounts 
cd $path/data/featureCounts
featureCounts -a $path/data/sequences/RNAs_100.saf -o $path/featureCounts_counts.txt 2> $path/featureCounts_info.txt -F "SAF" *.sam


