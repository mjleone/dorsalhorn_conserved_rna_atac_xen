#!/bin/bash

gwas_list=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/deepika_pain_gwas/pain_gwas_list.txt
while IFS= read -r line
do
  wc -l $line
done < $gwas_list
