#!/bin/bash
#SBATCH --job-name=SALA_Neuron_THP1_robust
#SBATCH --partition=batch
#SBATCH --output=gene0_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=24:00:00

perl /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/assemble_gene_annotator_v0.1.pl \
--chrom_size_path=/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/resources/v1.0.1/resources/genome/hg38.gencode_v39/tsv/chrom.sizes.tsv \
--model_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/Neuron_THP1.S3/bed/Neuron_THP1.S3.model.bed.bgz \
--model_info_gz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/Neuron_THP1.S3/log/Neuron_THP1.S3.model.info.tsv.gz \
--revert_ref_model_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.bed.bgz \
--ref_model_gene_link=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript_to_gene.tsv \
--novel_gene_prefix=IN1G \
--disable_ref_chain_bound_gene_anno=yes \
--min_ref_exon_overlap_pct=10 \
--exon_overlap_dist=-1 \
--locus_merge_dist=100000 \
--exclude_t_type=retained_intron \
--out_prefix=iPSC_NSC_Neuron.S3.t0.gene \
--out_dir=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/table0_gene \
--max_thread=1 \
--bedtools_bin=/usr/bin/bedtools \
--tabix_bin=/usr/bin/tabix \
--bgzip_bin=/usr/bin/bgzip

