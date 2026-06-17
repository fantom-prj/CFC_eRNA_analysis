#!/bin/bash
#SBATCH --job-name=gtf_Neuron_THP1_robust
#SBATCH --partition=batch
#SBATCH --output=gtf_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Path to the Perl script
SCRIPT="/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/code/others/SALA.gene_gtf_annotation.R"

# Activate base explicitly
source /home/yip/anaconda3/etc/profile.d/conda.sh
conda activate base

# Execute filter
/usr/bin/Rscript "$SCRIPT" \
/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/Neuron_THP1.S3_try \
Neuron_THP1.S3 \
/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/resources \
/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/resources/GENCODE_V39 \
/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/table4_gene/Neuron_THP1.S3


