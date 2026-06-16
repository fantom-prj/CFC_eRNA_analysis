#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed

# Read each line (sample name and paths) from the text file
while read -r sample fq_path; do
    echo "Processing sample: $sample"

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed/${sample}

    # Exclude supplementary alignment, Convert BAM to BED and gzip
    samtools view -F 2048 -b ${sample}_genome.sort.bam | bedtools bamtobed -bed12 -i stdin | gzip > ${sample}.bed12.bed.gz
    
    cd ..

done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/ont.fq.path.txt

