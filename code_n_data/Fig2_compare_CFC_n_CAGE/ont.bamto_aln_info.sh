#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed

# Read each line (sample name and paths) from the text file
while read -r sample path; do
    echo "Processing sample: $sample"

    # Remove supplementary and secondary alignments and Convert BAM to alignment info and gzip
    
    samtools view -F 256 -F 2048 "$path" | awk '{print $1, $5, gensub(/.*NM:i:([0-9]+).*/, "\\1", "g", $0)}' | gzip > "${sample}.tsv.gz"

    
done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed/ont.raw.bam.path.txt



samtools view input.bam | awk '{print $1, $5, gensub(/NM:i:(\d+)/, "\\1", "g", $0)}'