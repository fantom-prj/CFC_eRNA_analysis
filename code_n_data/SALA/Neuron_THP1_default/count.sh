#!/bin/bash
#SBATCH --job-name=count_Neuron_THP1_robust
#SBATCH --partition=batch
#SBATCH --output=SALA_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Path to the Perl script
COUNT_SCRIPT="/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/code/others/SALA.count_matrix.R"

# Execute the count
/usr/bin/Rscript "$COUNT_SCRIPT" \
/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/Neuron_THP1.S3 \
/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/Neuron_THP1_robust2/transcript/Neuron_THP1.S3/log \
/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/resources/GENCODE_V39/transcript_to_gene.tsv
