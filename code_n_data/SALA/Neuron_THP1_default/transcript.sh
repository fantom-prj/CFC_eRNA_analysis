#!/bin/bash
#SBATCH --job-name=SALA_Neuron_THP1_robust
#SBATCH --partition=batch
#SBATCH --output=transcript_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=24:00:00

# Path to the Perl script
SALA_SCRIPT="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/end5_guided_assembler_v0.1.20231102.pl"

# Execute the assembler
perl $SALA_SCRIPT \
--qry_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231110/Neuron_THP1/ontCAGE/bam_to_bed/pool_bed/Neuron_THP1.bed.bgz \
--ref_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.bed.bgz \
--chrom_size_path=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/tsv/chrom.sizes.tsv \
--out_dir=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript \
--max_thread=1 \
--out_prefix=Neuron_THP1.S3 \
--min_transcript_length=15 \
--doubtful_end_avoid_summit=yes \
--min_exon_length=1 \
--print_trnscrptID=no \
--chrom_fasta_path=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/fasta/genome.fa \
--min_output_qry_count=1 \
--trnscpt_set_end_priority=summit:commonest:longest \
--doubtful_end_merge_dist=150 \
--novel_model_prefix=ROBT \
--conf_end3_merge_flank=150 \
--conf_end5_merge_flank=75 \
--conf_end5_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1_S3/input_end3_end5_junct_bed/ontCAGE.Neuron_THP1.end5.cluster.bed.bgz \
--conf_end3_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/ontCAGE.Neuron_THP1.end3.cluster.bed.bgz \
--min_summit_dist_split=50 \
--retain_no_qry_ref_bound_set=no \
--min_size_split=100 \
--min_frac_split=0.2 \
--signal_end5_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1_S3/input_end3_end5_junct_bed/ontCAGE.Neuron_THP1.end5.signal.bed.bgz \
--signal_end3_bed_bgz=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/ontCAGE.Neuron_THP1.end3.signal.bed.bgz \
--conf_end3_add_ref=yes \
--conf_end5_add_ref=yes \
--conf_junction_bed=/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/Gencode_v39.junct.bed,/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/Neuron_THP1.long_read.hi_qual.junct.bed,/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/input_end3_end5_junct_bed/iPSC_NSC_Neuron.short_read.hi_qual.junct.bed \
--min_qry_score=0 \
--bedtools_bin=/usr/bin/bedtools \
--tabix_bin=/usr/bin/tabix \
--bgzip_bin=/usr/bin/bgzip

