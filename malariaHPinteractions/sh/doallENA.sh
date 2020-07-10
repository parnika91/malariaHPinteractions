#!/bin/bash

run=$1
echo $run

default_path="/SAN/Plasmo_compare/"

ulimit -m 130000000000
ulimit -v 130000000000

studyID="SRP171171"
host=mouse
parasite=Pyoelii

# ena links
bash --verbose $default_path/SRAdb/Scripts/malariaHPinteractions/sh/DownloadSRAtoolkit.sh $run # download run or pair of runs, based on layout=Paired or Single

# Step 2: map

#cd /SAN/Plasmo_compare/fastq_download_tmp/

#if [ -e $run\_1.fastq.gz ]; then
  # put .gz; then add --readFilesCommand zcat
#  STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$run_1.fastq.gz /SAN/Plasmo_compare/fastq_download_tmp/$run\_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$run --readFilesCommand zcat

#else
#  STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$run.fastq --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$run --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

#fi
#cd ..

#cp Aligned.out.bam $run.bam

# Step 3: gene expression quantification

##Rscript --vanilla $default_path/SRAdb/Scripts/malariaHPinteractions/R/countOverlaps.R $host $parasite $run $studyID --save

# Step 4: Enter used runs in blacklist

##cat $default_path/$studyID/runs.txt | grep $run >> $default_path/Input/blacklist.txt # also study
 # Step 5: remove fastq and bam files
##rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.fastq
#rm /SAN/Plasmo_compare/ncbi/sra/$run.sra
##rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.out
#rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.tab
#fi
##rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.bam



