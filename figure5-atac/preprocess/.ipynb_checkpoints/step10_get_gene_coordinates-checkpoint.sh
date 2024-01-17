#!/bin/bash

#conda activate base

PATH=$PATH:~/repos/gene-based-ldsc-functions

cd /projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/celltype_specific_enhancers/tables/gene-based-ldsc

bash map_genes.sh

cd ~/repos/gene-based-ldsc-functions/

GENE_DIR='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/celltype_specific_enhancers/tables/gene-based-ldsc'

for gene_set in $GENE_DIR/*csv
do
    echo $gene_set
   python get_gene_and_exon_coordinates.py -i "$gene_set"
    
done