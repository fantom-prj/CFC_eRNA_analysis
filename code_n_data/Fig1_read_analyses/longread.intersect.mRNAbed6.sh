#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed

# Read each line (sample name and paths) from the text file
while read -r sample fq_path; do
    echo "Processing sample: $sample"

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed/${sample}

    # bed12 intersect with mRNA bed6
    bedtools intersect -c -s -a ${sample}.bed12.bed.gz -b /analysisdata/fantom6/Interactome/resources/gencode.v39.annotation.mRNA.bed6.sort.bed -F 0.9 | gzip > ${sample}.mRNA.bed.gz
    cd ..

done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/ont.fq.path.txt

