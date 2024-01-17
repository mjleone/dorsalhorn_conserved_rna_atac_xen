#!/bin/bash


EnrichmentsDir='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/genebased_gwas_enrichments_mouse/'

path_to_fore_annots=${EnrichmentsDir}"annotations/DorsalHorn/foregroundAnnotations/"
path_to_back_annots=${EnrichmentsDir}"annotations/DorsalHorn/backgroundAnnotations/DHS_DorsalHorn_all_cell_types"

path_to_fore_bedfiles=${EnrichmentsDir}"foregroundBedFiles/*"

path_to_output=${EnrichmentsDir}"DorsalHorn_ldsc.ldcts"

# make file
> $path_to_output
for file in $path_to_fore_bedfiles
do 
	fore_name=$(basename ${file::-4})
	full_line="${fore_name}	${path_to_fore_annots}${fore_name}.,${path_to_back_annots}."
    # append to file
	echo $full_line >> $path_to_output
done
