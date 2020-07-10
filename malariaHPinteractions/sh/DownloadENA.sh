#!/bin/bash

fastq="$1"
echo $fastq

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh $fastq /SAN/Plasmo_compare/fastq_download_tmp/