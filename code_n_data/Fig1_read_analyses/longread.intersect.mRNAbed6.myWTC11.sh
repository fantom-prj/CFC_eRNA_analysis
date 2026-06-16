#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed/myWTC11

# Read each line (sample name and paths) from the text file
while read -r sample bed12_path; do
    echo "Processing sample: $sample"

    # bed12 intersect with mRNA bed6
    bedtools intersect -c -s -a ${bed12_path} -b /analysisdata/fantom6/Interactome/resources/gencode.v39.annotation.mRNA.bed6.sort.bed -F 0.9 | gzip > ${sample}.mRNA.bed.gz

done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/myWTC11.bed12.path.txt

