#!/bin/bash

# With this script, the URLs from download.txt are plugged in to download fastq files

fastq="$1"
echo $fastq

ascp -QT -l 300m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh $fastq tmp/