#!/bin/bash

# input files
studies='positive_experiments.txt'
getRuns='GetRuns.R'

# count number of studied = 3
number_studies=( `cat $studies | wc -l` )


for study in $(cat $studies); do
	
  studyID=$(echo $study | cut -d, -f 1)
  host=$(echo $study | cut -d, -f 2)
  parasite=$(echo $study | cut -d, -f 3)
  layout=$(echo $study | cut -d, -f 4)

  Rscript --vanilla $getRuns $studyID --save

  #runs='runs_rs.txt'
  cp runs_rs.txt runs_$studyID.txt

  # count number of lines in runs_rs.txt file
  #runNumber=$(cat runs_rs.txt | wc -l)

  # create folder for study
  mkdir $studyID
  mv runs_$studyID.txt -t $studyID
  cp $host$parasite.fa $host$parasite.fa.amb $host$parasite.fa.ann $host$parasite.fa.bwt $host$parasite.fa.pac $host$parasite.fa.sa -t $studyID # time consuming step
  #cp $parasite.fa.amb $parasite.fa.ann $parasite.fa.bwt $parasite.fa.pac $parasite.fa.sa -t $studyID
  #cp .fa $parasite.fa > $study
  cd $studyID

  #touch map_$studyID.txt
  #echo -e "Run\tHost_percent\tParasite_percent" >> map_$studyID.txt


  # if case for 0 Gb runs or no run IDs
  if [ ! -s runs_$studyID.txt ]; then
  echo "No runs available; leaving study and going to next study, if any"

  # fastq-dump for each run
  elif [ -f runs_$studyID.txt ]; then

    for line in $(cat runs_$studyID.txt);do 

      current=$line
      bash ../GetRuns$layout.sh $current 
      #echo -e "$line\t$host_per\t$par_per" >> map_$studyID.txt
      
      if [ -e $current\_1.fastq ]; then

        # bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
        bwa mem -t 15 $host$parasite.fa $current\_1.fastq $current\_2.fastq > aln_$host$parasite$current.sam
        samtools view -S -b aln_$host$parasite$current.sam > aln_$host$parasite$current.bam

      elif [ -e $current.fastq ]; then

        #bwa mem ref.fa reads.fq > aln-se.sam
        bwa mem -t 15 $host$parasite.fa $current.fastq > aln_$host$parasite$current.sam
        samtools view -S -b aln_$host$parasite$current.sam > aln_$host$parasite$current.bam

      fi

    done
  fi
  cd ..
done

