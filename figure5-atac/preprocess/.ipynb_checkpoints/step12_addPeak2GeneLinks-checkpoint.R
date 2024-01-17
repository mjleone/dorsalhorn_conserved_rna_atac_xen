# mike ran: conda activate base

suppressMessages(library(ArchR))
library(SingleCellExperiment)
library(SummarizedExperiment)
library(here); library(tidyverse)
library(Seurat)
library(rtracklayer)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = 8)
GENOME = 'mm10'

PLOTDIR='figures/exploratory/celltype_specific_enhancers/plots'
DATADIR='data/tidy_data/celltype_specific_enhancers'

#######################
## get ArchR project
PROJDIR=here("data/tidy_data/ArchRProjects/Mouse_DorsalHorn_scATAC")
proj = loadArchRProject(path = PROJDIR)

proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "IterativeLSI200",  # see step 1 of label_mouse
    corCutOff = 0.1
)


proj = saveArchRProject(ArchRProj = proj)
