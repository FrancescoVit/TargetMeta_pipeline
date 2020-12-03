#!/bin/bash

########################################################################################
########################################################################################

# This is a collection of bash command, used for the initial steps of data analysis;
# from raw sequences to OTUs. Those command are used for analysis of targeted 
# metagenomics data on the 16S

# Script author: Francesco Vitali
# Contact: francesco.vitali@ibba.cnr.it
# Github: FrancescoVit

# This pipeline assumes that the following programs are installed
# 
# 1. FASTQC (in path) https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# 2. MULTIQC (in path) https://multiqc.info/
# 3. CUTADAPT (not in path) https://cutadapt.readthedocs.io/en/stable/ 
# 4. SICKLE (in path) https://github.com/najoshi/sickle
# 5. MICCA (in path) https://micca.readthedocs.io/en/latest/ 
# 6. RDP classifier (not in path) https://github.com/rdpstaff/classifier 

# Credits for any of the above used programs goest to the respective authors. 

# DISCLAIMER: Use these command and the information contained here at your own risk!
# I'm not responsible for loss of data

# DISCLAIMER: this is a continuous work in progress, both on the side of tools udes, 
# and on the side of "scripting style". I have been working towards higher variable 
# use an parametrization. In some cases, absolute or relative file path are included
# and should be changed as necessary. Moreover, also some specific file names or 
# references to project name should be changed from time to time
# 
########################################################################################
########################################################################################


## Set some parameters for the script
currpath=$(pwd) # top project directory
core=8 # number of core to use

## Create analysis folder
mkdir $currpath/Raw_reads
mkdir $currpath/raw_quality 
mkdir $currpath/renamed/
mkdir $currpath/cutadapt/
mkdir $currpath/quality_control/
mkdir $currpath/quality_control/cutadapt
mkdir $currpath/quality_control/trimmed
mkdir $currpath/trimmed
mkdir $currpath/MICCA_16S_WP2

## Visualize quality of raw reads
fastqc $currpath/Raw_reads*.fastq.gz -o ./raw_quality/ -t 8
cd $currpath/raw_quality 
multiqc .
cd $currpath

## Rename samples. We usually have files from a service provider, with the usual not useful information
## as L001 or S1. Here cutting to the actual sample name and read number

cd $currpath/Raw_reads

este=".fastq.gz"

for i in *.fastq.gz
do
  sample=$(echo "$i" | cut -d "_" -f1)
  read=$(echo "$i" | cut -d "_" -f4)
  cp $i $currpath/renamed/$sample"_"$read$este
  echo -e "$i\t-->\t$sample"_"$read$este" >> log_renamer.txt
done

cd ..

## Count reads

cd $currpath/renamed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc " " >> seq_count_16S_raw.txt
done

## Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2
for i in *.fastq.gz
do 
echo "$i" | cut -d "_" -f1 >> names.txt
sed 'n; d' names.txt > names_single.txt  
done

cp names_single $currpath/trimmed

# remove primers with cutadapt
# Primers for V3-V4 16S (Sequences from Pindo and FMACH)
# forward: CCTACGGGNGGCWGCAG
# reverse: GACTACNVGGGTWTCTAATCC

cd $currpath/renamed

while read file
do
	echo "Running cutadapt on file "${i}""
	/usr/local/bin/cutadapt -g Forward=CCTACGGGNGGCWGCAG -G Reverse=GACTACNVGGGTWTCTAATCC --discard-untrimmed --pair-filter=any -o $currpath/cutadapt/"${file}_R1_cutadapt.fastq.gz" -p $currpath/cutadapt/"${file}_R2_cutadapt.fastq.gz" "${file}_R1.fastq.gz" "${file}_R2.fastq.gz" >> $currpath/quality_control/cutadapt/cutadapt_report.txt  
done < names_single.txt

# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any

cd $currpath/cutadapt

fastqc ./*.fastq.gz -o $currpath/quality_control/cutadapt -t 8

cd  $currpath/quality_control/cutadapt
multiqc .


# All R2 reads have a drop in quality around 280 bp, consider clipping befor trimming
cd $currpath/cutadapt

while read file
do
	echo "Running cutadapt on R2 of "${file}""
	/usr/local/bin/cutadapt -l 275 -o "${file}_R2_short_cutadapt.fastq.gz" "${file}_R2_cutadapt.fastq.gz" >> $currpath/quality_control/cutadapt/cutadapt_report_shorteningR2.txt
done < names_single.txt

# trim low quality part in 5'

cd $currpath/cutadapt

while read file
do
	echo "Running sickle on file "${file}""
	echo "Running sickle on file "${file}"" >> $currpath/quality_control/trimmed/stats_trim.txt
	sickle pe -f "${file}_R1_cutadapt.fastq.gz" -r "${file}_R2_short_cutadapt.fastq.gz" -o $currpath/trimmed/"${file}_trimmed_R1.fastq.gz" -p $currpath/trimmed/"${file}_trimmed_R2.fastq.gz" -s $currpath/trimmed/"${file}_singles.gz" -t sanger -q 25 -g 1>> $currpath/quality_control/trimmed/stats_trim.txt
done < names_single.txt

cd $currpath/trimmed
fastqc ./*.fastq.gz -o $currpath/quality_control/trimmed -t 8


cd  $currpath/quality_control/trimmed
multiqc .

# Multiqc says: all green!

## Count reads in the raw data and in the final QC data to highlight loss NOT WORKING

#echo "Name" "R1_raw" "R1_QC" "R2_raw" "R2_QC" > seq_count_16S_allQC.txt

# print header

#while read i
#do
#        echo "Counting "${i}""
#        echo -n $i >> seq_count_16S_allQC.txt
#        echo -n ""
#        echo -n $(zcat $currpath/renamed/"${i}_R1.fastq.gz" | wc -l) / 4 | bc >> seq_count_16S_allQC.txt
#        echo -n ""
#        echo -n $(zcat $currpath/trimmed/"${i}_trimmed_R1.fastq.gz" | wc -l) / 4 | bc >> seq_count_16S_allQC.txt
#        echo -n ""
#        echo -n $(zcat $currpath/renamed/"${i}_R2.fastq.gz" | wc -l) / 4 | bc >> seq_count_16S_allQC.txt
#        echo -n ""
#        echo $(zcat $currpath/trimmed/"${i}_trimmed_R2.fastq.gz" | wc -l) / 4 | bc >> seq_count_16S_allQC.txt
#done < names_single.txt


# Count QC treated data

cd $currpath/trimmed

for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc " " >> seq_count_16S_QC.txt
done

# Read loss from raw to QC is: 3.6% mean, 6.3% max and 3.05% min

## Join reads (MICCA): https://micca.readthedocs.io/en/latest/

cd $currpath/trimmed

# remove singles reads from sickle

rm ./*singles.gz

gunzip *.fastq.gz

micca mergepairs -i $currpath/trimmed/*_R1.fastq -o $currpath/MICCA_16S_WP2/WP2_assembled_16S.fastq -l 100 -d 8 -t 7

# -l : minimum overlap between reads 
# -d : maximum mismatch in overlap region

# Counting reads in assembled file

grep -c '^@M' $currpath/MICCA_16S_WP2/WP2_assembled_16S.fastq   

# 5277248. Read loss from QC reads to assembled is  17%


# Count assembled reads per sample

while read i
do
	echo " " >> seq_count_assembled.txt
	echo -n $i " " >> seq_count_assembled.txt
	grep -wc "$i" WP2_assembled_16S.fastq >> seq_count_assembled.txt
done < names_single.txt


# Remove N from assembly
micca filter -i $currpath/MICCA_16S_WP2/WP2_assembled_16S.fastq -o $currpath/MICCA_16S_WP2/WP2_assembled_16S.fasta --maxns 0

# Count
grep -c '>' $currpath/MICCA_16S_WP2/WP2_assembled_16S.fasta


# to improve picking, also for UNOISE maybe, but surely for DADA2, reads should be of same length. From FASTQC on the assebled, a good minimal length would be 400 bp

/usr/local/bin/cutadapt -l 400 -m 400 -o $currpath/MICCA_16S_WP2/WP2_assembled_16S_400bp.fasta $currpath/MICCA_16S_WP2/WP2_assembled_16S.fasta


# pick otu
micca otu -m denovo_unoise -i $currpath/MICCA_16S_WP2/WP2_assembled_16S_400bp.fasta -o $currpath/MICCA_16S_WP2/ -t 8 --rmchim

# classify RDP

export RDPPATH=/home/fvitali/Documenti/Personal_PATH_folder/rdp_classifier_2.13/
micca classify -m rdp -i $currpath/MICCA_16S_WP2/otus.fasta --rdp-gene 16srrna -o $currpath/MICCA_16S_WP2/taxa.txt


# plot stats, using the file before N removal, hence there are some little differences. The command needs a fastq, but after N removal I have a fasta 

micca stats -i $currpath/MICCA_burkina_transmic/ITS_burkina_assembled.fastq -o stats_full



##################################################################################################################################################################