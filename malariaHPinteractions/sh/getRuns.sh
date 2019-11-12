#!/bin/bash

default_path="/SAN/Plasmo_compare/SRAdb"
#RgetRuns="$default_path/Scripts/GetRuns.R"

#declare -a current_studies

current_studies=$1

cd /SAN/Plasmo_compare/SRAdb/

Rscript $default_path/Scripts/malariaHPinteractions/R/GetRuns.R $current_studies

#for i in ${current_studies[@]}; do
#  Rscript --vanilla $RgetRuns $i --save
#  mkdir $default_path/Output/$i
#  mv runs_rs.txt runs_$i.txt # leave runs of this study outside of the study folder so that the downloads can happen easily
#  cp runs_$i.txt -t $default_path/Output/$i
  
#done
# rm runs_*.txt
# rm /SAN/Plasmo_compare/SRAdb/Input/current_runs.txt

#touch Input/current_runs.txt
#for i in ${current_studies[@]}; do # runs_rs.txt has run and study IDs
#  cat $default_path/Output/$i/runs_$i.txt >> $default_path/Input/current_runs.txt
#done


