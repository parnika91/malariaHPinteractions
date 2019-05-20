#!/bin/bash

hp=$1

default_path="/SAN/Plasmo_compare/SRAdb/"
studies='positive_experiments.txt'

current_studies=()

for study in $(cat $default_path/Input/$studies | cut -f 1); do
	
  studyID=$study
  host=$(grep $study $default_path/Input/$studies | cut -f 2)
  parasite=$(grep $study $default_path/Input/$studies | cut -f 3)
  #layout=$(echo $study | cut -d, -f 4)

  if [ $host$parasite = $hp ]; then
    current_studies=("$studyID" "${current_studies[@]}")
    #current_layout=("$layout" "${current_layout[@]}")
  fi

done
echo "I am looking for $hp"

#mkdir /SAN/Plasmo_compare/Genomes/indices/$hp

# check if hp pair has been indexed

cd /SAN/Plasmo_compare/Genomes/indices/$hp

if [ -f genomeParameters.txt ]; then
  echo "Indexed files exist"
else
  echo "Indexing $hp.fa now..."
    #index="STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /SAN/Plasmo_compare/Genomes/indices/$hp --genomeFastaFiles /SAN/Plasmo_compare/Genomes/fasta/{1} --sjdbGTFfile /SAN/Plasmo_compare/Genomes/annotation/{2} --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --limitGenomeGenerateRAM 210000000000" # 195GB=210000000000bytes
    # index
    # Got gtf files, you dont need to specify gene parent and transcript parent
    index="STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /SAN/Plasmo_compare/Genomes/indices/$hp --genomeFastaFiles /SAN/Plasmo_compare/Genomes/fasta/{1} --sjdbGTFfile /SAN/Plasmo_compare/Genomes/annotation/{2} --limitGenomeGenerateRAM 210000000000"
    parallel --xapply $index {1} {2} ::: $hp.fa ::: $hp.gtf
fi

#STAR --runThreadN 12 --genomeDir /SAN/Plasmo_compare/SRAdb/Genomes/fasta/ --genomeFastaFiles /SAN/Plasmo_compare/SRAdb/Genomes/fasta/mousePyoelii.fa --sjdbGTFfile /SAN/Plasmo_compare/SRAdb/Genomes/annotation/mousePyoelii.gff3 --sjdbOverhang 100 --readFilesIn /SAN/Plasmo_compare/fastq.download.tmp/*.fastq --sjdbGTFtagExonParentTranscript Parent --sjdbGTFtagExonParentGene Parent --twopassMode Basic

# STAR --runMode genomeGenerate --genomeDir /SAN/Plasmo_compare/SRAdb/Genomes/indices/mousePyoelii --genomeFastaFiles /SAN/Plasmo_compare/SRAdb/Genomes/fasta/mousePyoelii.fa --runThreadN 12

cd $default_path
