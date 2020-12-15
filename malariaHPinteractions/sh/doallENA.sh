#!/bin/bash

# With this script, for every run, the reads are mapped and sent to an R script for counting

run_study=$1
echo $run_study
runID="$(echo $run | cut -d, -f1)"
studyID="$(echo $run | cut -d, -f2)"

# take studyID, host and parasite info from positive_experiments.txt file.

positive_experiments="positive_experiments.txt"

host="$(grep $studyID $positive_experiments | cut -d, -f2)"
parasite="$(grep $studyID $positive_experiments | cut -d, -f3)"

sh getHPstudies.sh $host$parasite

cd tmp/

echo "Time to map run $runID of study $studyID!"

if [ -e $runID\_1.fastq.gz ]; then # mapping paired end reads
    STAR --runThreadN 6 \
    --genomeDir Genomes/indices/$host$parasite \
    --readFilesIn tmp/$runID\_1.fastq.gz tmp/$runID\_2.fastq.gz \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix tmp/$runID \
    --readFilesCommand zcat

  else # mapping single end reads
    STAR --runThreadN 6 \
    --genomeDir Genomes/indices/$host$parasite \
    --readFilesIn tmp/$runID.fastq \
    --outFileNamePrefix tmp/$runID \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat
fi

# moving with info about unique mapping percentage, etc and splice junctions into the study folder

cat tmp/$runID*.final.out > tmp/$runID\_$studyID.final.out
mv tmp/$runID\_$studyID.final.out $studyID/
mv tmp/$runID*.out.tab -t $studyID/

# Gene expression quantification

Rscript --vanilla countOverlaps.R $host $parasite $runID $studyID --save

# Enter used runs in "finished-runs" list

if [ -e $studyID/countWithGFF3_$runID.txt ]; then
  cat $studyID/runs\_$studyID.txt | grep $run_study >> finished_runs.txt # also study
 
  # Removing fastq and bam files to save space

  rm tmp/$runID*.fastq.gz
  rm tmp/$runID*.bam

fi