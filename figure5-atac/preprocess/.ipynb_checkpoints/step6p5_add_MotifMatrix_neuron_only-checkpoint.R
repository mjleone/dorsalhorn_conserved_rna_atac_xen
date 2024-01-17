
## Run below in terminal
#inter_pfen1 50 12
#condact archr_r4
#Rscript --no-save --no-restore --verbose step6c_add_MotifMatrix_neuron_only.R > logs/step6c.Rout 2>&1

library(ArchR)
library(parallel)
library(rhdf5)
#library(SingleCellExperiment)
#library(SeuratDisk)
#library(Seurat)
library(cowplot)
library(here)
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

##################################
### set Arrow File parameters ####
#addArchRThreads(threads = round(parallel::detectCores()*1/4))
addArchRThreads(threads = 8)

##################################
### load rheMac10 ArchR genome ###
# GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

addArchRGenome("mm10")


PROJDIR=here('data/tidy_data/ArchRProjects')
ARCHDIR=file.path(PROJDIR,'Mouse_scATAC_DorsalHorn_neuron2_wDeviations')
proj = loadArchRProject(ARCHDIR)

library(BSgenome.Mmusculus.UCSC.mm10)

########## needs to be archr and R version 4!

#BiocManager::install("JASPAR2020")
library(JASPAR2020)

proj <- addDeviationsMatrix(proj,  peakAnnotation = "Motif", force = TRUE)
proj = saveArchRProject(proj)
print('done!')