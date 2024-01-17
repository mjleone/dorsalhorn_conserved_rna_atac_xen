#!/bin/bash  

# install atac data pipeline: 

SEARCH_DIR='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/raw_data/peak_mm10'
#echo $SEARCH_DIR

COUNTER=0

for narrowPeak in $SEARCH_DIR/*
do
    #echo "$narrowPeak"
    
    #if (($COUNTER > 9))
    #then
    
    #echo "$COUNTER"
    
    sbatch -w compute-1-11 --mem 2000 ~/repos/atac_data_pipeline/scripts/halper_map_peak_orthologs.sh \
    -s Mus_musculus \
    -t Homo_sapiens \
    -o /scratch/mleone2/Mouse_DH_toHomo_sapiens/ \
    -b $narrowPeak \
    --keepChrPrefix chr \
    --halPath /projects/pfenninggroup/install/halLiftover
    
    #fi
    
    
    COUNTER=$((COUNTER+1))
    #echo "$COUNTER"
done