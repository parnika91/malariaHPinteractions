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

url=$prefix/$dir1/$dir2/$accession/


# download fastq files in /tmp/

# cd /tmp

if [ ! -f $fastq.fastq ]; then

  wget -c $url/$fastq*.fastq.gz -P /tmp
  gunzip /tmp/$fastq*.fastq.gz
fi

