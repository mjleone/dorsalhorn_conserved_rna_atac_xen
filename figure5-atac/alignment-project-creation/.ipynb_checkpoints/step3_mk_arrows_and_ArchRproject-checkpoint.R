suppressMessages(library(ArchR))
library(here)
library(tidyverse)

ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

###############################################
### get the fragments file from snATAC-seq ####
fragments_fn = here('data/raw_data/bed') %>%
  list.files(pattern = '.bed.gz$', full.names = T)
names(fragments_fn) = fragments_fn %>% basename() %>% ss('\\.', 1)

########################################
## load genomeAnnotation, geneAnnotation
addArchRThreads(threads = round(parallel::detectCores()/2))
GENOMEDIR='/home/bnphan/resources/genomes/rheMac10'
# contains `geneAnnotation` and `genomeAnnotation` objects
load(file.path(GENOMEDIR,'rheMac10_liftoff_GRCh38.p13_ArchR_annotations.rda'))

#############################################################################
### make that arrow file from bgzipped fragments file (tsv.gz or bed.gz) ####
ArrowFiles <- createArrowFiles( inputFiles = fragments_fn, 
                                sampleNames = names(fragments_fn), 
                                minTSS = 4, #Dont set this too high because you can always increase later
                                minFrags = 1000,
                                addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                geneAnnotation = geneAnnotation, #must be the custom rheMac10 version
                                genomeAnnotation = genomeAnnotation, #must be the custom rheMac10 version
                                force = TRUE)

doubleScores <- addDoubletScores(input = ArrowFiles, 
                                 k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                                 LSIMethod = 1,  knnMethod = "UMAP") #Refers to the embedding to use for nearest neighbor search with doublet projection.

###############################
### make the ArchR Project ####
proj = ArchRProject( ArrowFiles = ArrowFiles,
                     outputDirectory = here("data/tidy_data/ArchRProjects/Mouse_DorsalHorn_scATAC"), copyArrows = TRUE,                geneAnnotation = geneAnnotation, #must be the custom rheMac10 version
                     genomeAnnotation = genomeAnnotation #must be the custom rheMac10 version
)

## filter doublets
proj = filterDoublets( proj, cutEnrich = .5, cutScore = -Inf, filterRatio = 1)

## increase unique fragments cutoff to 10^3.5, remove cluster of low QC cell 
idxSample <- BiocGenerics::which(proj$nFrags > 10^3.5)
cellsSample <- proj$cellNames[idxSample]
proj = subsetCells(ArchRProj = proj, cellNames = cellsSample)

## increase TSSEnrichment cutoff to 8, remove cluster of low QC cell 
idxSample2 <- BiocGenerics::which(proj$TSSEnrichment > 8 )
cellsSample2 <- proj$cellNames[idxSample2]
proj = subsetCells(ArchRProj = proj, cellNames = cellsSample2)

## Add metadata columns
proj$Tissue = 'dorsalHorn'

table(proj$Sample) 
proj = saveArchRProject(proj)

############################################
### move Arrow files and QC to data dir ####
rsync = 'rsync -Paq --remove-source-files'
files = c(paste0(names(fragments_fn),'.arrow'), 'QualityControl', 'ArchRLogs')
dir = here("data/raw_data/arrow")
thecall = paste(rsync, files, dir)
sapply(thecall, system)

