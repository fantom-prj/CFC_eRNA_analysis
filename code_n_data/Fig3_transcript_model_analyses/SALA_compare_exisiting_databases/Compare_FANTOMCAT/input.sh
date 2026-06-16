#!/bin/bash
cwd=`dirname $0`
cd $cwd

transcript_bed_to_end_bed_bigwig='/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/assembly/20231114/test_overlay_FCAT/input/00_transcript_bed_to_end_bed_bigwig.pl'
assemble_table4_bed_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/all_gtf_file/table5pENST.bed12.bed.gz'
assemble_table1_model_info_tsv_gz='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/table4_gene/Neuron_THP1_T4_10percent/log/Neuron_THP1_T4_10percent.model.info.tsv.gz'
FANTOM_CAT_permissive_transcript_bed_gz='/home/hon-chun/resources/genome/human/inUse/hg38/FANTOM_CAT/CAT_liftoverer/00_output/lv2_permissive.hg38/bed/lv2_permissive.trnscpt.hg38.bed.gz'
FANTOM_CAT_permissive_DPIClstr_bed_gz='/home/hon-chun/resources/genome/human/inUse/hg38/FANTOM_CAT/CAT_liftoverer/00_output/lv2_permissive.hg38/bed/lv2_permissive.DPIClstr.hg38.bed.gz'
chrom_size='/home/hon-chun/resources/genome/human/inUse/hg38/fasta/hg38.chrom.sizes.sorted.txt'
slop=150
tag='Neuron_THP1_FCAT'

zcat $assemble_table1_model_info_tsv_gz | tail -n +2 | awk -F'\t' '{print $1"\t"$7"\t"$9"\t"$9"\t"$8}' >./$tag.table1.transcript_to_gene.tsv
zcat $assemble_table4_bed_gz | bgzip -c >./$tag.table4.bed.bgz
tabix -p bed ./$tag.table4.bed.bgz

zcat $FANTOM_CAT_permissive_transcript_bed_gz | bgzip -c >./FANTOM_CAT_permissive_transcript.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive_transcript.bed.bgz

zcat $FANTOM_CAT_permissive_DPIClstr_bed_gz | bgzip -c >./FANTOM_CAT_permissive_DPIClstr.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive_DPIClstr.bed.bgz

$transcript_bed_to_end_bed_bigwig ./FANTOM_CAT_permissive_transcript.bed.bgz $chrom_size $slop FANTOM_CAT_permissive $cwd

cat ./FANTOM_CAT_permissive.end3.flank.with_summit.bed | bgzip -c >./FANTOM_CAT_permissive.end3.flank.with_summit.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive.end3.flank.with_summit.bed.bgz

cat ./FANTOM_CAT_permissive.end5.flank.with_summit.bed | bgzip -c >./FANTOM_CAT_permissive.end5.flank.with_summit.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive.end5.flank.with_summit.bed.bgz

cat ./FANTOM_CAT_permissive.end3.bed | bgzip -c >./FANTOM_CAT_permissive.end3.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive.end3.bed.bgz

cat ./FANTOM_CAT_permissive.end5.bed | bgzip -c >./FANTOM_CAT_permissive.end5.bed.bgz
tabix -p bed ./FANTOM_CAT_permissive.end5.bed.bgz
