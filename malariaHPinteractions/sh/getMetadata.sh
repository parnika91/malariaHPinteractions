#!/bin/bash
positive_experiments="/SAN/Plasmo_compare/SRAdb/Input/positive_experiments.txt"

for study in $(cat $positive_experiments); do
  studyID=$(echo $study | cut -d',' -f1)
  wget -O ./$studyID_info.csv 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= $studyID' > /SAN/Plasmo_compare/SRAdb/Output/$study/$studyID_info.csv
done
