#!/bin/bash
#SBATCH -p pfen1
#SBATCH --mem=80G
#SBATCH --error=logs/subcluster_%A_output.txt
#SBATCH --output=logs/subcluster_%A_output.txt


source ~/.bashrc
#conda activate r4
conda activate ArchR_Env

Rscript step3_subcluster_neuron_glia.R
