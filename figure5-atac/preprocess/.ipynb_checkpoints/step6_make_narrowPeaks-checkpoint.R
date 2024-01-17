ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
suppressMessages(library(rtracklayer))
library(tidyverse)
library(here)

source('/projects/pfenninggroup/machineLearningForComputationalBiology/snATAC_cross_species_caudate/code/raw_code/hal_scripts/narrowPeakFunctions.R')

##############################
### read in ArchR project ####
DATADIR='data/tidy_data'
CODEDIR='code/raw_code/label_mouse_snATAC_cells'

##############################
### read in ArchR project ####
PROJDIR='data/tidy_data'
LABEL='Mouse_DH'; GENOME = 'mm10'; 

#######################################################
# get the reproducible peaks across Cluster2 cell types
peak_rds_fn = list.files(path = here(DATADIR,'ArchRProjects', 'Mouse_DorsalHorn_scATAC','PeakCalls'), 
                           full.names = T, pattern = '.rds')
names(peak_rds_fn) = ss(basename(peak_rds_fn),'-reproduciblePeaks.gr.rds')
# only want the reproducible peaks
peak_rds_fn <- mylist[grep('-reproduciblePeaks.gr.rds', peak_rds_fn)]

peakList = lapply(peak_rds_fn, readRDS)



# label the summit and peak name using rheMac10 coordinates
peakList = lapply(peakList, addSummitCenter)
peakList = lapply(peakList, nameNarrowPeakRanges, genome = GENOME)
peakList = lapply(peakList, sort)

###############################################
# create directory and narrowPeak file names ##
PEAKDIR2=here('data/raw_data','peak_mm10')
system(paste('mkdir -p',PEAKDIR2))
narrowPeak_mm10_fn = here(PEAKDIR2, paste0(LABEL, '.', names(peakList), 
                                          '.mm10.narrowPeak.gz'))
# write peaks to narrowPeak file
outRanges = mapply(write_GRangesToNarrowPeak,
                   gr = peakList, file = narrowPeak_mm10_fn, genome = GENOME)

