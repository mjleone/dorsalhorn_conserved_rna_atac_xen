#### delete the content of PathToAnnotOutput before re-running
#### I commented out the first part. bedtools is ok on head node?
#### I have my own version of annotate_bed_LDSC

#!/bin/bash

conda activate ldsc

#PathToNarrowPeak="/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/raw_data/halper/*narrowPeak.gz"
PathToBed="/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/celltype_specific_enhancers/tables/gene-based-ldsc/introns_and_flanks/*.bed.gz"
PathToBackgroundBed="/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/celltype_specific_enhancers/tables/gene-based-ldsc-background/introns_and_flanks/*.bed.gz"

EnrichmentsDir='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/genebased_gwas_enrichments_mouse/'

if [ ! -d "$EnrichmentsDir" ]; then
  mkdir "$EnrichmentsDir"
  echo "made $EnrichmentsDir"
fi

BackgroundDir=${EnrichmentsDir}"backgroundBedFiles"
if [ ! -d "$BackgroundDir" ]; then
  mkdir "$BackgroundDir"
  echo "made $BackgroundDir"
fi

PathToOutput=${EnrichmentsDir}"backgroundBedFiles/DorsalHorn_all_cell_types.bed"

echo "Making union of foregrounds"

### make the union of all foregrounds - should take a few seconds
zcat ${PathBed} | awk 'OFS="\t" {print $1, $2, $3}'| bedtools sort -i stdin | bedtools merge -i stdin > ${PathToOutput}
gzip ${PathToOutput}

echo "Done. Merging union of foregrounds with DHS"

### merge the union of all foregrounds with DHS - should take a few seconds

PathToMerged=${EnrichmentsDir}"backgroundBedFiles/DorsalHorn_all_cell_types.bed.gz"
#PathToDHS="/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/mouseDHS/mm10_univ_dhs_ucsc_flank300_halLifted2hg38_orthologs.bed.gz"

# Use BackgroundBed (genes) instead of DHS
PathToOutput=${EnrichmentsDir}"backgroundBedFiles/DHS_DorsalHorn_all_cell_types.bed"
zcat ${PathToMerged} ${PathToBackgroundBed} | awk 'OFS="\t" {print $1, $2, $3}' | bedtools sort -i stdin | bedtools merge -i stdin > ${PathToOutput}


echo "Done. Annotating background"

### annotate the background - should take a few hours (around 3)
AnnotationsTopFolder=${EnrichmentsDir}"annotations/"
if [ ! -d "$AnnotationsTopFolder" ]; then
  mkdir "$AnnotationsTopFolder"
  echo "made $AnnotationsTopFolder"
fi


Annotations=${EnrichmentsDir}"annotations/DorsalHorn/"
if [ ! -d "$Annotations" ]; then
  mkdir "$Annotations"
  echo "made $Annotations"
fi

PathToAnnotOutput=${EnrichmentsDir}"annotations/DorsalHorn/backgroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${PathToOutput} -n DHS_DorsalHorn_all_cell_types -o ${PathToAnnotOutput}

echo "Done. Converting to bed"
### convert narrowPeak foreground files to bed file format - should take a few seconds

PathToOutput=${EnrichmentsDir}"foregroundBedFiles/"
if [ ! -d "$PathToOutput" ]; then
  mkdir "$PathToOutput"
  echo "made $PathToOutput"
fi

cp ${PathToBed} ${PathToOutput}
gunzip ${PathToOutput}*

#for file in ${PathToNarrowPeak}; do x=${file##*/}; cat ${PathToOutput}${x/.narrowPeak.gz/.bed}; done
#for file in ${PathToBed}; do x=${file##*/}; cat ${file} | awk 'OFS="\t" {print $1, $2, $3}' > ${PathToOutput}${x}; done

#echo "Done. Annotating foregrounds"

#### annotate foregrounds - should take a few hours (around 6+)
PathToAnnotOutput=${EnrichmentsDir}"annotations/DorsalHorn/foregroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

#for file in ${PathToOutput}*.bed; do echo ${file##*/}; done

#for file in ${PathToOutput}*.bed;
for file in ${PathToOutput}*.bed; do x=${file##*/}; base_file_without_bed=${x/.bed}; sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i $file -n $base_file_without_bed -o ${PathToAnnotOutput}; done

# munge gwas'
# create ldcts file
# run: sbatch file_name.sb