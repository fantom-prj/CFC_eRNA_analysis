#!/bin/bash
#SBATCH --job-name=ABC_3celltypes
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=48:00:00
#SBATCH --output=logs/ABC_%j.out
#SBATCH --error=logs/ABC_%j.err

set -eo pipefail

cd /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main
mkdir -p logs

source /home/yip/anaconda3/etc/profile.d/conda.sh
conda activate abc-env

snakemake \
  --cores ${SLURM_CPUS_PER_TASK} \
  --rerun-incomplete \
  --keep-going \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/iPSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/NSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  /home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/Neuron/Predictions/EnhancerPredictionsAllPutative.tsv.gz


