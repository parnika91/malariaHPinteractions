#!/bin/bash

for studyID in SRP189198; do
#studyID=$1
#runs="/SAN/Plasmo_compare/SRAdb/Output/$studyID/runs\_$studyID.txt"

#current_runs="/SAN/Plasmo_compare/SRAdb/Input/current_runs.txt"

#cut -d, -f 1 current_runs.txt | parallel --eta -j 12 --link bash $doall {1}

##awk 'FS="\t", OFS="\t" { gsub("ftp.sra.ebi.ac.uk", "era-fasp@fasp.sra.ebi.ac.uk:"); print }' /SAN/Plasmo_compare/SRAdb/Output/$studyID/$studyID\_accession.txt | cut -f7 | awk -F ";" 'OFS=" " {print $1, $2}' | awk NF > /SAN/Plasmo_compare/SRAdb/Output/$studyID/download.txt

# ls /SAN/Plasmo_compare/fastq_download_tmp/ | grep fastq | rev | cut -d '_' -f3- | rev | uniq > /SAN/Plasmo_compare/Kai/macrophage_runs.txt
 #parallel --eta -j 3 --link bash --verbose /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/sh/doall.sh ::: $(cut -d, -f 1 /SAN/Plasmo_compare/SRAdb/Output/$studyID/download.txt)
  
parallel --eta -j 16 --link bash --verbose /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/sh/DownloadENA.sh ::: $(cat /SAN/Plasmo_compare/SRAdb/Output/$studyID/download.txt)

##wait

##parallel --eta -j 2 --link bash --verbose /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/sh/doallENA.sh ::: $(cat /SAN/Plasmo_compare/SRAdb/Output/$studyID/download.txt)
  parallel --eta -j 3 --link bash --verbose /SAN/Plasmo_compare/SRAdb/Scripts/malariaHPinteractions/sh/doallENA.sh ::: $(cat /SAN/Plasmo_compare/SRAdb/Output/$studyID/runs\_$studyID.txt)
done

exit
