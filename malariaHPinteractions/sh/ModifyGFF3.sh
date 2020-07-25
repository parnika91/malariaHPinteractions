#! /bin/bash/


# remove all comment lines
grep -v '^#.' Pknowlesi.gtf > temp && mv temp Pknowlesi.gtf

# attach prefix to gtf
sed -e 's/^/Pkn_chr/' Pknowlesi.gtf > temp && mv temp Pknowlesi.gtf


# attach prefix to fasta
sed -e 's/>/>Pkn_chr/' Pknowlesi.fa > temp && mv temp Pknowlesi.fa

