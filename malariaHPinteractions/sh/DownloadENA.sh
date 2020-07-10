#!/bin/bash
# from http://www.nxn.se/valent/streaming-rna-seq-data-from-ena
fastq="$1"
echo $fastq

# adding ascp prefix to each line of download.txt

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh $fastq /SAN/Plasmo_compare/fastq_download_tmp/

#/localstorage/parnika/sratoolkit.2.10.8-ubuntu64/bin/prefetch -a "/home/parnika/.aspera/connect/bin/ascp|/home/parnika/.aspera/connect/etc/asperaweb_id_dsa.openssh" $run
#sleep $(( RANDOM/300 ))
#/localstorage/parnika/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump -I --split-3 $run -O /SAN/Plasmo_compare/fastq_download_tmp/
