#!/bin/bash
cwd=`dirname $0`
cd $cwd

base_tag='Neuron_THP1'
ref_tag='refseq2024'

transcript_bed_to_end_bed_bigwig='./00_transcript_bed_to_end_bed_bigwig.pl'
assemble_bed_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/all_gtf_file/table5pENST.bed12.bed.gz'
assemble_info_tsv_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/table4_gene/Neuron_THP1_T4_10percent/log/Neuron_THP1_T4_10percent.model.info.tsv.gz'
ref_bed_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Compare_refseq2024/GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.gz'
base_rng_end5_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/Neuron_THP1.S3.end5.bed.bgz'
base_rng_end3_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/Neuron_THP1.S3.end3.bed.bgz'

merge_d_end5=75
merge_d_end3=150
base_rng_end5=./$base_tag.end5.bed
base_rng_end3=./$base_tag.end3.bed

zcat $base_rng_end5_gz | cut -f 1-6 >$base_rng_end5
zcat $base_rng_end3_gz | cut -f 1-6 >$base_rng_end3

zcat $assemble_info_tsv_gz | tail -n +2 | awk -F'\t' '{print $1"\t"$7"\t"$9"\t"$9"\t"$8}' >./$base_tag.transcript_to_gene.tsv
zcat $assemble_bed_gz | bgzip -c >./$base_tag.bed.bgz
tabix -p bed ./$base_tag.bed.bgz

zcat $ref_bed_gz  | sort -k1,1 -k2,2n | bgzip -c >./$ref_tag.bed.bgz
tabix -p bed ./$ref_tag.bed.bgz

zcat ./$ref_tag.bed.bgz ./$base_tag.bed.bgz | sort -k1,1 -k2,2n | bgzip -c >./$ref_tag.$base_tag.transcript.bed.bgz
tabix -p bed ./$ref_tag.$base_tag.transcript.bed.bgz

#$transcript_bed, $merge_d_end5, $merge_d_end3, $base_rng_end5, $base_rng_end3, $outPrefix, $outDir
perl $transcript_bed_to_end_bed_bigwig ./$ref_tag.$base_tag.transcript.bed.bgz $merge_d_end5 $merge_d_end3 $base_rng_end5 $base_rng_end3 $ref_tag.$base_tag $cwd

