#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed

# Read each line (sample name and paths) from the text file
while read -r sample path; do
    echo "Processing sample: $sample"

    # Remove supplementary reads and Convert BAM to BED and gzip
    
    samtools view -F 2048 -b "$path" | bedtools bamtobed -i stdin | gzip > ${sample}.bed.gz
    
done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed/ont.raw.bam.path.txt
