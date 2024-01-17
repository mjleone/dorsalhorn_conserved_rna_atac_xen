#### delete the content of PathToAnnotOutput before re-running
#### I commented out the first part. bedtools is ok on head node?
#### I have my own version of annotate_bed_LDSC

#!/bin/bash

conda activate ldsc

#PathToNarrowPeak="/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/raw_data/halper/*narrowPeak.gz"
PathToNarrowPeak='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/pfenliver_halper/MouseLiver_idr.optimal_peak.Mus_musculusToHomo_sapiens.HALPER.narrowPeak.gz'
EnrichmentsDir='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/gwas_enrichments_mouse_pfenliver/'

if [ ! -d "$EnrichmentsDir" ]; then
  mkdir "$EnrichmentsDir"
  echo "made $EnrichmentsDir"
fi

BackgroundDir=${EnrichmentsDir}"backgroundBedFiles"
if [ ! -d "$BackgroundDir" ]; then
  mkdir "$BackgroundDir"
  echo "made $BackgroundDir"
fi

PathToOutput=${EnrichmentsDir}"backgroundBedFiles/all.bed.gz"

### make the union of all foregrounds - should take a few seconds
cp $PathToNarrowPeak $PathToOutput


### merge the union of all foregrounds with DHS - should take a few seconds

PathToMerged=${EnrichmentsDir}"backgroundBedFiles/all.bed.gz"
PathToDHS="/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/background/mouseDHS/mm10_univ_dhs_ucsc_flank300_halLifted2hg38_orthologs.bed.gz"

PathToOutput=${EnrichmentsDir}"backgroundBedFiles/DHS_all.bed"
zcat ${PathToMerged} ${PathToDHS} | awk 'OFS="\t" {print $1, $2, $3}' | bedtools sort -i stdin | bedtools merge -i stdin > ${PathToOutput}

### annotate the background - should take a few hours (around 3)

Annotations=${EnrichmentsDir}"annotations/"
if [ ! -d "$Annotations" ]; then
  mkdir "$Annotations"
  echo "made $Annotations"
fi

PathToAnnotOutput=${EnrichmentsDir}"annotations/backgroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i ${PathToOutput} -n DHS_all -o ${PathToAnnotOutput}

### convert narrowPeak foreground files to bed file format - should take a few seconds

PathToOutput=${EnrichmentsDir}"foregroundBedFiles/"
if [ ! -d "$PathToOutput" ]; then
  mkdir "$PathToOutput"
  echo "made $PathToOutput"
fi

# make empty files; not, overwrites
#for file in ${PathToNarrowPeak}; do x=${file##*/}; cat ${PathToOutput}${x/.narrowPeak.gz/.bed}; done
for file in ${PathToNarrowPeak}; do x=${file##*/}; zcat ${file} | awk 'OFS="\t" {print $1, $2, $3}' > ${PathToOutput}${x/.narrowPeak.gz/.bed}; done

### annotate foregrounds - should take a few hours (around 6+)
PathToAnnotOutput=${EnrichmentsDir}"annotations/foregroundAnnotations/"
if [ ! -d "$PathToAnnotOutput" ]; then
  mkdir "$PathToAnnotOutput"
  echo "made $PathToAnnotOutput"
fi

#for file in ${PathToOutput}*.bed; do echo ${file##*/}; done

for file in ${PathToOutput}*.bed; do x=${file##*/}; base_file_without_bed=${x/.bed}; sbatch --partition pool1 /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/scripts/annotate_bed_LDSC_hg38.sh -i $file -n $base_file_without_bed -o ${PathToAnnotOutput}; done

# munge gwas'
# create ldcts file
# run: sbatch file_name.sb