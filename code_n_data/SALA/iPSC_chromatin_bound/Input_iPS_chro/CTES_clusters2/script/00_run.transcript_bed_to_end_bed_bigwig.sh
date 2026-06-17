#!/bin/bash
cwd=`dirname $0`
cd $cwd
outDir='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTES_clusters2/end3_bed_bigwig'
outPrefix='iPSchro'
transcript_bed_to_end_bed_bigwig='./transcript_bed_to_end_bed_bigwig.pl'
transcript_bed='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231110/iPSchro/ontCAGE/bam_to_bed/pool_bed_all/iPSchro.bed.bgz'
chrom_size_path='/home/hon-chun/resources/genome/human/inUse/hg38/fasta/hg38.chrom.sizes.sorted.txt'
perl $transcript_bed_to_end_bed_bigwig $transcript_bed $chrom_size_path $outPrefix $outDir
cp $0 $outDir