#!/bin/sh

base_dir="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed"

# Read each line (sample name and paths) from the text file
while read -r sample r1_path r2_path; do
    echo "Processing sample: $sample"

    mkdir -p "$base_dir/$sample"
    mkdir -p "$base_dir/$sample/PE_mapping"
    cd "$base_dir/$sample/PE_mapping"

    # Run STAR for paired-end reads
    /home/scratch/moirai/bin/STAR \
        --twopassMode None \
        --genomeDir /home/scratch/moirai/STAR/hg38analysisset \
        --readFilesIn "$r1_path" "$r2_path" \
        --outSAMmultNmax -1 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --outReadsUnmapped Fastx \
        --runThreadN 16 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000

    # Convert BAM to BED and gzip
    samtools view -b -f 0x1 -F 0x100 Aligned.sortedByCoord.out.bam | bedtools bamtobed -i stdin | gzip > "${sample}.primary.PE.R1.bed.gz"

	mkdir -p "$base_dir/$sample/SE_mapping"
	cd "$base_dir/$sample/SE_mapping"

    # Run STAR for single-end reads
    /home/scratch/moirai/bin/STAR \
        --twopassMode None \
        --genomeDir /home/scratch/moirai/STAR/hg38analysisset \
        --readFilesIn "$r1_path" \
        --outSAMmultNmax -1 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --outReadsUnmapped Fastx \
        --runThreadN 16 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000

    # Convert BAM to BED and gzip
    samtools view -b -F 0x100 Aligned.sortedByCoord.out.bam | bedtools bamtobed -i stdin | gzip > "${sample}.primary.SE.R1.bed.gz"

    cd "$base_dir"  # Go back to the base directory

done < "$base_dir/fq.path.txt"


#--outSAMmultNmax -1 is default
#extract all the NH:i: tag from the bam. this represent the total loci being mapped of a read.
#samtools view -h /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/repeat_element/bamtobed/iPS_rep1/Aligned.sortedByCoord.out.bam | awk -F '\t' '{for(i=12; i<=NF; ++i) if($i ~ /^NH:/) print $i}' > ssips1.nh_tags.txt
