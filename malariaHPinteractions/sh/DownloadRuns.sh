#!/bin/bash
# from http://www.nxn.se/valent/streaming-rna-seq-data-from-ena
fastq="$1"
echo $fastq

prefix=ftp://ftp.sra.ebi.ac.uk/vol1/fastq

accession=$(echo $fastq | tr '.' '_' | cut -d'_' -f 1)

dir1=${accession:0:6}

a_len=${#accession}
if (( $a_len == 9 )); then
    dir2="";
elif (( $a_len == 10 )); then
    dir2=00${accession:9:1};
elif (( $a_len == 11)); then
    dir2=0${accession:9:2};
else
    dir2=${accession:9:3};
fi

url=$prefix/$dir1/$dir2/$accession/

# download fastq files in /tmp/

cd /SAN/Plasmo_compare/fastq_download_tmp

if [ ! -f $fastq.fastq.gz ]; then
  
  content=$(w3m -dump $url | grep '_1')

  if [ ! -z "$content" ];
  then
      wget -c $url/$fastq\_1.fastq.gz
      gunzip $fastq\_1.fastq.gz
      wget -c $url/$fastq\_2.fastq.gz
      gunzip $fastq\_2.fastq.gz

      if [ ! -f $fastq\_1.fastq ]; then
        /SAN/Plasmo_compare/SRAdb/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch -a "/home/parnika/.aspera/connect/bin/ascp|/home/parnika/.aspera/connect/etc/asperaweb_id_dsa.openssh" $fastq
        /localstorage/parnika/PapersAndBooks/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump -I --split-files $fastq
      fi
  else
      wget -c $url/$fastq.fastq.gz
      gunzip $fastq.fastq.gz
      if [ ! -f $fastq.fastq ]; then
        /SAN/Plasmo_compare/SRAdb/sratoolkit.2.8.2-1-ubuntu64/bin/prefetch -a "/home/parnika/.aspera/connect/bin/ascp|/home/parnika/.aspera/connect/etc/asperaweb_id_dsa.openssh" $fastq
        /localstorage/parnika/PapersAndBooks/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump -I --split-files $fastq
      fi
  fi
  #wget -c -A gz -r -l 1 -nd $url
fi
