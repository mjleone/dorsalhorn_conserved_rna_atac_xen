#!/bin/bash

PathToAnnotations='/projects/pfenninggroup/bloodAD/ldsc/IDR_optimal_peaks_hg38/annotations/'
EnrichmentsDir='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq/data/tidy_data/gwas_enrichments_immunecalvin/'

if [ ! -d "$EnrichmentsDir" ]; then
  mkdir "$EnrichmentsDir"
  echo "made $EnrichmentsDir"
fi

Annotations=${EnrichmentsDir}"annotations/"
if [ ! -d "$Annotations" ]; then
  mkdir "$Annotations"
  echo "made $Annotations"
fi

PathToBackAnnotOutput=${EnrichmentsDir}"annotations/backgroundAnnotationsRoadmapAll/"
if [ ! -d "$PathToBackAnnotOutput" ]; then
  mkdir "$PathToBackAnnotOutput"
  echo "made $PathToBackAnnotOutput"
fi

PathToForeAnnotOutput=${EnrichmentsDir}"annotations/foregroundAnnotations/"
if [ ! -d "$PathToForeAnnotOutput" ]; then
  mkdir "$PathToForeAnnotOutput"
  echo "made $PathToForeAnnotOutput"
fi

backAnnotations="$PathToAnnotations"Roadmap_All*
cp $backAnnotations $PathToBackAnnotOutput

#leave out roadmap and archive folder
foreAnnotations="$PathToAnnotations"*.*
cp $foreAnnotations $PathToForeAnnotOutput
rm "$PathToForeAnnotOutput"annot.sh
rm "$PathToForeAnnotOutput"*Roadmap*
rm "$PathToForeAnnotOutput"count.ipynb




#####
##### copy the bed files
#####
####
foreBed=${EnrichmentsDir}"foregroundBedFiles/"
if [ ! -d "$foreBed" ]; then
  mkdir "$foreBed"
  echo "made $foreBed"
fi

backBed=${EnrichmentsDir}"backgroundBedFilesRoadmapAll/"
if [ ! -d "$backBed" ]; then
  mkdir "$backBed"
  echo "made $backBed"
fi

peaksDir="/projects/pfenninggroup/bloodAD/ldsc/IDR_optimal_peaks_hg38/All_peaks/"
roadmap_all_bed="/projects/pfenninggroup/bloodAD/ldsc/IDR_optimal_peaks_hg38/All_peaks/Roadmap_All_merged.bed"

cp $roadmap_all_bed "$backBed"Roadmap_All_merged.bed
cp "$peaksDir"* $foreBed

rm "$foreBed"cat_all.ipynb
rm "$foreBed"*Roadmap*

gunzip "$foreBed"*.narrowPeak.gz
mv "$foreBed"*.narrowPeak "$foreBed"*.narrowPeak

for file in $foreBed*.narrowPeak
do 
    mv $file "${file::-11}".bed
done

#### NOTE: I had to manually change the THP1Snyder_macrophages bed file name to match the annotations. It was previously THP1Snyder_macrophageSnyder



