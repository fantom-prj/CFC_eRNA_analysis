#!/bin/bash
cwd=`dirname $0`
cd /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Compare_refseq2024/sala

ref_tag='refseq2024'
end5_guided_assembler_script='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/end5_guided_assembler_v0.1.updated.pl'

perl $end5_guided_assembler_script \
--qry_bed_bgz=./input/$ref_tag.bed.bgz \
--ref_bed_bgz=./input/Neuron_THP1.bed.bgz \
--conf_end5_bed_bgz=./input/$ref_tag.Neuron_THP1.end5.merge.filter.bed.bgz \
--conf_end3_bed_bgz=./input/$ref_tag.Neuron_THP1.end3.merge.filter.bed.bgz \
--signal_end5_bed_bgz=./input/$ref_tag.Neuron_THP1.end5.signal.bed.bgz \
--signal_end3_bed_bgz=./input/$ref_tag.Neuron_THP1.end3.signal.bed.bgz \
--out_prefix=$ref_tag.Neuron_THP1 \
--out_dir=./transcript \
--max_thread=1 \
--min_transcript_length=25 \
--doubtful_end_avoid_summit=yes \
--min_exon_length=1 \
--print_trnscrptID=no \
--conf_junction_bed=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/Gencode_v39.junct.bed,/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/Neuron_THP1.long_read.hi_qual.junct.bed,/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/iPSC_NSC_Neuron.short_read.hi_qual.junct.bed \
--chrom_size_path=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/tsv/chrom.sizes.tsv \
--chrom_fasta_path=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/fasta/genome.fa \
--min_output_qry_count=1 \
--trnscpt_set_end_priority=commonest:summit:longest \
--doubtful_end_merge_dist=150 \
--conf_end3_merge_flank=150 \
--conf_end5_merge_flank=75 \
--min_summit_dist_split=50 \
--novel_model_prefix=NEWT \
--print_trnscrptID=yes \
--retain_no_qry_ref_bound_set=yes \
--doubtful_end_avoid_summit=yes \
--min_size_split=100 \
--min_frac_split=0.2 \
--conf_end3_add_ref=yes \
--conf_end5_add_ref=yes \
--min_qry_score=0 \
--bedtools_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bedtools/bedtools \
--tabix_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/tabix/tabix \
--bgzip_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bgzip/bgzip

