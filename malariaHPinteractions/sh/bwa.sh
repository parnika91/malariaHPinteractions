#! /bin/bash

# to perform counts

# get the number of folders with SRA experiments
declare -a folders
folders=( `find -type d -name '*RP*' -print` )

# get rid of the './' for the folder names
for fn in ${folders[@]}; do ${fn:2}; done

# check if there are bam files in these folders, and if there are any, use Rscript to count
for folder in "${folders[@]}"; do
  if [ ! -e *.bam ]; then

    host=$(ls | cat '*.fa.amb' | cut -dP -f 1)
    parasite=P$(ls | echo *.fa.amb | cut -dP -f 2 | cut -d. -f 1)

    Rscript --vanilla ../countOverlaps.R $host $parasite --save

  else
    echo "Nothing to map!"
  fi
