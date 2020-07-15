#!/bin/bash
run=$1
echo $run

default_path="/SAN/Plasmo_compare/"


studyID="ERP109432"
host=human
parasite=Pberghei

#example: era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR826/004/SRR8263174/SRR8263174_1.fastq.gz
runID="$(cut -d, -f1 <<< $run)"

ulimit -m 130000000000
ulimit -v 130000000000

# Step 2: map

cd /SAN/Plasmo_compare/fastq_download_tmp/

if [ -e $runID\_1.fastq.gz ]; then
#put .gz; then add --readFilesCommand zcat
  STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID\_1.fastq.gz /SAN/Plasmo_compare/fastq_download_tmp/$runID\_2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --readFilesCommand zcat

else
  STAR --runThreadN 6 --genomeDir /SAN/Plasmo_compare/Genomes/indices/$host$parasite --readFilesIn /SAN/Plasmo_compare/fastq_download_tmp/$runID.fastq --outFileNamePrefix /SAN/Plasmo_compare/fastq_download_tmp/$runID --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

fi
cd ..

# Step 3: gene expression quantification

Rscript --vanilla $default_path/SRAdb/Scripts/malariaHPinteractions/R/countOverlaps.R $host $parasite $runID $studyID --save

# Step 4: Enter used runs in blacklist

cat $default_path/$studyID/runs.txt | grep $run >> $default_path/Input/blacklist.txt # also study
# Step 5: remove fastq and bam files
rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.fastq
rm /SAN/Plasmo_compare/fastq_download_tmp/$run*.bam



