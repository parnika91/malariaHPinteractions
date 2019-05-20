#!/bin/bash

# rename Plasmodium and host .fa and .gff3 files so that the .bam files have same names as the .gff3 files while counting

declare -a fasta_files
fasta_files=$(ls *\_ens.fa)

for file in ${fasta_files[@]}; do sed -i 's/>/>'${file:0:3}'_chr/g' $file; done

grep 'CHR_HSCHR' human_test.fa | sed -i 's/>CHR_HSCHR/>Hs_chr'${file:0:3}'/g'

for i in $(grep 'CHR_HSCHR' human_test.fa); sed -i 's/>CHR_HSCHR/>Hs_chr'${i:9:1}'/g'; done

for i in $(grep 'CHR_HSCHR' human_test.fa); do
  name=$(echo $i | cut -d ' ' -f 1)
  sed -i 's/>CHR_HSCHR/>Hs_chr'${name:9:1}'/g' human_test.fa;
done

