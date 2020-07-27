#!/bin/bash

run=$1
echo $run
runID="$(echo $run | cut -d, -f1)"
studyID="$(echo $run | cut -d, -f2)"

# take studyID

default_path="/SAN/Plasmo_compare/SRAdb"
positive_experiments="/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt"

host="$(grep $studyID $positive_experiments | cut -d, -f2)"
parasite="$(grep $studyID $positive_experiments | cut -d, -f3)"

cd /SAN/Plasmo_compare/fastq_download_tmp/


if [ $host$parasite == "monkeyPknowlesi" ]; then
  if [ -e $runID\_1.fastq.gz ]; then
    /SAN/Plasmo_compare/SRAdb/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID\_1.fastq.gz /SAN/Plasmo_compare/fastq_download_tmp/$runID\_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --readFilesCommand zcat

  else
    /SAN/Plasmo_compare/SRAdb/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID.fastq --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
  fi
else
  if [ -e $runID\_1.fastq.gz ]; then
    /SAN/Plasmo_compare/SRAdb/STAR-2.6.0c/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID\_1.fastq.gz /SAN/Plasmo_compare/fastq_download_tmp/$runID\_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --readFilesCommand zcat

  else
    /SAN/Plasmo_compare/SRAdb/STAR-2.6.0c/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID.fastq --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat
  fi
fi

cat /SAN/Plasmo_compare/fastq_download_tmp/$runID*.final.out > /SAN/Plasmo_compare/fastq_download_tmp/$runID\_$studyID.final.out
mv /SAN/Plasmo_compare/fastq_download_tmp/$runID\_$studyID.final.out /SAN/Plasmo_compare/SRAdb/Output/$studyID/
mv /SAN/Plasmo_compare/fastq_download_tmp/$runID*.out.tab -t /SAN/Plasmo_compare/SRAdb/Output/$studyID/

# Step 3: gene expression quantification

Rscript --vanilla /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/R/countOverlaps.R $host $parasite $runID $studyID --save

# Step 4: Enter used runs in blacklist

if [ -e $default_path/Output/$studyID/countWithGFF3_$runID.txt ]; then
  cat /SAN/Plasmo_compare/SRAdb/Output/$studyID/runs\_$studyID.txt | grep $run >> /SAN/Plasmo_compare/SRAdb/Input/blacklist.txt # also study
 
# Step 5: remove fastq and bam files
  rm /SAN/Plasmo_compare/fastq_download_tmp/$runID*.fastq.gz
  rm /SAN/Plasmo_compare/fastq_download_tmp/$runID*.bam

fi