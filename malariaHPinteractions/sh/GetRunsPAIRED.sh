#!/bin/bash
# from http://www.nxn.se/valent/streaming-rna-seq-data-from-ena
fastq="$1"

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

url1=$prefix/$dir1/$dir2/$accession/$fastq\_1.fastq.gz
url2=$prefix/$dir1/$dir2/$accession/$fastq\_2.fastq.gz

if [ ! -f $fastq.fastq ]; then

  wget -c $url1 -P /tmp
  wget -c $url2 -P /tmp
  gunzip /tmp/$fastq\_1.fastq.gz
  gunzip /tmp/$fastq\_2.fastq.gz
fi
