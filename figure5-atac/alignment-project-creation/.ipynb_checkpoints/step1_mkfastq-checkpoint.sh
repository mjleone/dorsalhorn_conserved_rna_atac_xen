#!/bin/bash
#SBATCH --partition=pool1
#SBATCH --time=1-0:00:00
#SBATCH --mem=20G
#SBATCH --export=ALL
#SBATCH --signal=2
#SBATCH --job-name=download
#SBATCH --error=logs/download_%A_%a_out.txt
#SBATCH --output=logs/download_%A_%a_out.txt
#SBATCH --no-requeue

TMPDIR=/scratch/bnphan
PROJDIR=/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq
CODEDIR=$PROJDIR/code/raw_code/preprocess_mouse_snATACseq
DATADIR=$PROJDIR/data/raw_data

mkdir -p $DATADIR/fastq; cd $DATADIR/

cellranger-atac mkfastq \
--id=Mouse_SealDorsalHorn_snATAC \
--localcores=4 --localmem=15 --run=$DATADIR/bcl \
--samplesheet=${DATADIR}/tables/SampleSheet.csv \
--output-dir=$DATADIR/fastq

# conserved sequence signature improve association