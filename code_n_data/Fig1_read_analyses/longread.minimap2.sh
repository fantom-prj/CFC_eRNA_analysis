#!/bin/sh

cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed

# Read each line (sample name and paths) from the text file
while read -r sample fq_path; do
    echo "Processing sample: $sample"

mkdir ${sample}
cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/bamtobed/${sample}

    # Run minimap2 for the current sample
    ~/minimap2-2.26_x64-linux/minimap2 -t 10 -2 -ax splice -uf -k 15 \
    --MD \
    --junc-bed /osc-fs_home/scratch/callum/Genome_files/gencode.v39.annotation.bed \
    --secondary=no /osc-fs_home/scratch/callum/Genome_files/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.new.mmi \
    ${fq_path} | samtools sort -o ${sample}_genome.sort.bam --write-index - 

    # Exclude supplementary alignment, Convert BAM to BED and gzip
    samtools view -F 2048 -b ${sample}_genome.sort.bam | bedtools bamtobed -i stdin | gzip > ${sample}.bed.gz
    
    cd ..

done < /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/other.full.length.long.read/ont.fq.path.txt

