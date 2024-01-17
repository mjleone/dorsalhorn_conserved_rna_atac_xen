#!/bin/bash
 
#remember to activate hal

narrowPeak="/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/pfenliver_halper/MouseLiver_idr.optimal_peak.narrowPeak"

sbatch -w compute-1-11 --mem 2000 ~/repos/atac_data_pipeline/scripts/halper_map_peak_orthologs.sh \
    -s Mus_musculus \
    -t Homo_sapiens \
    -o /scratch/mleone2/pfenliver_MusMusculus_toHomo_sapiens/ \
    -b $narrowPeak \
    --keepChrPrefix chr \
    --halPath /projects/pfenninggroup/install/halLiftover
   

