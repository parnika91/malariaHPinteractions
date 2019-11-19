#!/bin/bash

# get the input run from outside

run=$1
default_path="/SAN/Plasmo_compare/"

#if grep -q $run $default_path/Input/blacklist.txt; then exit; fi

#studies=$default_path/Input/positive_experiments.txt
ulimit -m 130000000000
ulimit -v 130000000000
##
##
##
##
# Step 1: download fastq file with runID

studyID=Macrophage_Kai
host=human
parasite=Pberghei


# bash --verbose $default_path/Scripts/malariaHPinteractions/sh/DownloadSRAtoolkit.sh $run # download run or pair of runs, based on layout=Paired or Single

# Step 2: map

cd /SAN/Plasmo_compare/fastq_download_tmp/

if [ -e $run\_R1_001.fastq ]; then
  # put .gz; then add --readFilesCommand zcat
  STAR --runThreadN 4 --genomeDir /SAN/Plasmo_compare/Genomes/indices/humanPberghei --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$run\_R1_001.fastq /SAN/Plasmo_compare/fastq_download_tmp/$run\_R1_001.fastq --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$run

else
  STAR --runThreadN 4 --genomeDir /SAN/Plasmo_compare/Genomes/indices/humanPberghei --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$run.fastq --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$run --outSAMtype BAM SortedByCoordinate

fi
cd ..

#cp Aligned.out.bam $run.bam

#if [ -e $run\_1.fastq ]; then
 # bwa mem -t 6 ../Genomes/fasta/mousePyoelii.fa /SAN/Plasmo_compare/fastq_download_tmp/$run\_1.fastq /SAN/Plasmo_compare/fastq_download_tmp/$run\_2.fastq > $run.aln.sam
#else
 # bwa mem -t 6 ../Genomes/fasta/mousePyoelii.fa /SAN/Plasmo_compare/fastq_download_tmp/$run.fastq > $run.aln.sam
#fi

#samtools view -S -b $run.aln.sam > $run.bam
#cat /SAN/Plasmo_compare/fastq_download_tmp/$run*.final.out > /SAN/Plasmo_compare/fastq_download_tmp/$run\_$studyID.final.out
mv /SAN/Plasmo_compare/fastq_download_tmp/*.final.out /SAN/Plasmo_compare/Kai/Mapping/
mv /SAN/Plasmo_compare/fastq_download_tmp/$run*.out.tab -t /SAN/Plasmo_compare/Kai/Mapping/

# Step 3: gene expression quantification

Rscript --vanilla $default_path/SRAdb/Scripts/malariaHPinteractions/R/countOverlaps.R $host $parasite $run $studyID --save

# Step 4: Enter used runs in blacklist

#if [ -e $default_path/Output/$studyID/countWithGFF3_$run.txt ]; then
 #cat $default_path/Input/current_runs.txt | grep $run >> $default_path/Input/blacklist.txt # also study
 # Step 5: remove fastq and bam files
rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.fastq
#rm /SAN/Plasmo_compare/ncbi/sra/$run.sra
rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.out
#rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.tab
#fi
rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.bam



