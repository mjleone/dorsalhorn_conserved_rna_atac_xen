#!/bin/bash
#SBATCH -p pfen1
#SBATCH --mem=80G
#SBATCH --error=logs/step1_output.txt
#SBATCH --output=logs/step1_output.txt


source ~/.bashrc
conda activate r4 
Rscript step1_harmony_clusters.R

