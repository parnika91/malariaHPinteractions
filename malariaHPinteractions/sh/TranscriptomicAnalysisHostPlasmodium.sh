#!/bin/bash

# declare -a hosts
# declare -a parasites

# hosts=("human" "mouse" "monkey")
# parasites=("Pberghei" "Pvivax" "Pfalciparum" "Pchabaudi" "Pyoelii" "Pcoatneyi" "Pcynomolgi")

# input files

host=human
parasite=Pberghei

default_path="/SAN/Plasmo_compare/"
#getHPstudies=$default_path/Scripts/getHPstudies.sh
getRuns=$default_path/Scripts/malariaHPinteractions/sh/getRuns.sh
doall=$default_path/Scripts/malariaHPinteractions/sh/doall.sh

#for i in $(cat $default_path/Input/positive_experiments.txt); do
 # host=$(echo $i | cut -d, -f 2)
 # parasite=$(echo $i | cut -d, -f 3)
 # echo $host$parasite >> $default_path/Input/hp.txt
#done

#cat $default_path/Input/hp.txt | sort | uniq > $default_path/Input/hp.txt.tmp && mv $default_path/Input/hp.txt.tmp $default_path/Input/hp.txt

#echo "$( parallel --gnu echo {1}{2} :::: hosts.txt :::: parasites.txt )\n" > hp.txt
#head -c -3 hp.txt > hp.txt.tmp && mv hp.txt.tmp hp.txt


# check if a study belongs to a host-parasite pair. If yes, download its runs. Else, move on to next study

# while loop begins here to do the whole thing again for another hp pair

#for pair in $(cat $default_path/Input/hp.txt); do

  #source $getHPstudies $pair
  #tail -n +2 hp.txt > hp.txt.tmp && mv hp.txt.tmp hp.txt # remove first hp pair


  #source $getRuns SRP108356 #${current_studies[@]}

  # cut -d, -f 1 current_runs.txt | parallel --eta -j 12 --link bash $doall {1}

  ls /SAN/Plasmo_compare/fastq_download_tmp/ | grep fastq | rev | cut -d '_' -f3- | rev | uniq > /SAN/Plasmo_compare/Kai/macrophage_runs.txt
  parallel --eta -j 3 --link bash --verbose /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/sh/doall.sh ::: $(cut -d, -f 1 /SAN/Plasmo_compare/Kai/macrophage_runs.txt)

  # Step 6: Final count tables and analysis

  #for i in ${current_studies[@]}; do
  # Rscript --vanilla $default_path/SRAdb/Scripts/malariaHPinteractions/R/count_analysis.R $host $parasite $i
  #  cp Distribution_$i.pdf -t $default_path/Output/$i
  #  cp Histogram_host_$i.pdf -t $default_path/Output/$i
  #  cp Histogram_parasite_$i.pdf -t $default_path/Output/$i
   ## cp table_$i.txt -t $default_path/Output/$i
   # cp hp_percent_$i.txt -t $default_path/Output/$i
  #done
  

# tempAnalysis: Rscript --vanilla tempAnalysis.R $host $parasite $study
#done

exit
