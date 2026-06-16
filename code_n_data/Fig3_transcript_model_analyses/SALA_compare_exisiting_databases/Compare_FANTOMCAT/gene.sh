#!/bin/bash
perl /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/assemble_gene_annotator_v0.1.pl \
--chrom_size_path=/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/resources/v1.0.1/resources/genome/hg38.gencode_v39/tsv/chrom.sizes.tsv \
--model_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/add_FANTOMCAT/sala/transcript/Neuron_THP1_FCAT/bed/Neuron_THP1_FCAT.model.bed.bgz \
--model_info_gz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/add_FANTOMCAT/sala/transcript/Neuron_THP1_FCAT/log/Neuron_THP1_FCAT.model.info.tsv.gz \
--revert_ref_model_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/all_gtf_file/table5pENST.bed12.bed.bgz \
--ref_model_gene_link=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/add_FANTOMCAT/sala/gene/table5pENST.info.tsv \
--out_prefix=Neuron_THP1_FCAT_disable_yes_10percent \
--out_dir=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/add_FANTOMCAT/sala/gene/ \
--novel_gene_prefix=NEWG \
--disable_ref_chain_bound_gene_anno=yes \
--min_ref_exon_overlap_pct=10 \
--exon_overlap_dist=-1 \
--locus_merge_dist=100000 \
--exclude_t_type=retained_intron \
--max_thread=13 \
--bedtools_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bedtools/bedtools \
--tabix_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/tabix/tabix \
--bgzip_bin=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bgzip/bgzip