#!/bin/bash
cwd=`dirname $0`
cd $cwd
outDir='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig'
outPrefix='Neuron_THP1'
transcript_bed_to_end_bed_bigwig='./transcript_bed_to_end_bed_bigwig.pl'
transcript_bed='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231110/Neuron_THP1/ontCAGE/bam_to_bed/pool_bed/Neuron_THP1.bed.bgz'
chrom_size_path='/home/hon-chun/resources/genome/human/inUse/hg38/fasta/hg38.chrom.sizes.sorted.txt'
perl $transcript_bed_to_end_bed_bigwig $transcript_bed $chrom_size_path $outPrefix $outDir
cp $0 $outDir