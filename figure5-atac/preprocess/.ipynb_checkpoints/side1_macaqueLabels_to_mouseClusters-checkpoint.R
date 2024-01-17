library(Seurat)
library(SeuratDisk)
library(harmony)
library(SingleCellExperiment)
library(tidyverse)
library(here)


##### This script separates the glia and neuron in the macaque snRNA dataset and also converts to mouse gene names, in preparation for ATAC-RNA integration

#### must run before step4a

###########################################################
## 0) # Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = x , mart = human, attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  
  # get 1-1 gene orthologs by ensembl
  oneToOne = split(genesV2$HGNC.symbol, genesV2$MGI.symbol) %>% 
    lengths()
  oneToOne = names(oneToOne)[which(oneToOne ==1)]
  genesV2 = genesV2[genesV2$MGI.symbol %in% oneToOne,]
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}

#######################################################
## 1) load in the neuron subtypemacaque Seurat objects  
obj_rm = LoadH5Seurat(here('data/tidy_data/rdas/macaque_integrated_neuronsv3.h5seurat'), misc = F, tools = F)
head(obj_rm[[]])
table(obj_rm$cluster_type) ## only contains neurons

## see if counts is data
tmp = all(obj_rm@assays$RNA@counts == obj_rm@assays$RNA@data) #yup...

## create new count object changing human gene names to 1-1 mouse ortholog
rm_to_mm_genes <- convertHumanGeneList(rownames(obj_rm))
counts_rna = round(2^obj_rm@assays$RNA@counts)
counts_rna = counts_rna[rownames(counts_rna) %in% rm_to_mm_genes$HGNC.symbol,]
rownames(counts_rna) = rm_to_mm_genes$MGI.symbol[match(rownames(counts_rna), 
                                                        rm_to_mm_genes$HGNC.symbol)]
obj_rm[["RNA"]] <- CreateAssayObject(counts = counts_rna)

## save this version of the dataset
SaveH5Seurat(obj_rm, misc = F, tools = F, overwrite = T,
             filename = here('data/tidy_data/rdas/macaque_integrated_neuronsv3.mouseGenes.h5seurat'))

sce_rm = as.SingleCellExperiment(obj_rm)
saveRDS(sce_rm, here('data/tidy_data/rdas/macaque_integrated_neuronsv3.mouseGenes.sce.rds'))



#######################################################
## 2) load in the neuron subtypemacaque Seurat objects  
sce_rm = readRDS(here('data/tidy_data/rdas/full_integrated_macaque_dorsal_hornv2.sce.rds'))
head(colData(sce_rm))

sce_rm$cell_type = make.names(as.character(sce_rm$cell_type))
sce_rm = sce_rm[,!grepl('GLUT|GABA|midVen', sce_rm$cell_type)] # drop the neurons
table(sce_rm$cell_type) ## no more neurons

## see if counts is same as logcounts
tmp = all(assays(sce_rm)$counts == assays(sce_rm)$logcounts) # not the same

## create new count object changing human gene names to 1-1 mouse ortholog
rm_to_mm_genes <- convertHumanGeneList(rownames(sce_rm))
counts_rna = assays(sce_rm)$counts
counts_rna = counts_rna[rownames(counts_rna) %in% rm_to_mm_genes$HGNC.symbol,]
rownames(counts_rna) = rm_to_mm_genes$MGI.symbol[match(rownames(counts_rna), 
                                                       rm_to_mm_genes$HGNC.symbol)]

sce_rm2 <- SingleCellExperiment(
  assays = list(counts = counts_rna), colData=colData(sce_rm)
)

## save this version of the dataset
saveRDS(sce_rm2, here('data/tidy_data/rdas/glia_integrated_macaque_dorsal_hornv2.mouseGenes.sce.rds'))


#######################################################
## 3) load in the subset of mouse snRNA-seq that Andreas made w/ human gene names 
obj_adultMouse = LoadH5Seurat(here('data/tidy_data/mmData.comp.adult.2fixerr.h5seurat'))

## find the raw mouse snRNA-seq data w/ mouse genes
h5ad_file = "/projects/pfenninggroup/singleCell/MacaqueMouse_SealDorsalHorn_snRNA-Seq/mmData/final_cluster_assignment.raw.mm.h5ad"
h5Seurat_file = here('data/tidy_data/rdas/final_cluster_assignment.raw.mm.h5Seurat')
Convert(h5ad_file, dest = h5Seurat_file, assay = "RNA",  overwrite = FALSE,  verbose = TRUE)
obj_mmData <- LoadH5Seurat(h5Seurat_file, misc = F, tools = F)

## subset to the same cells as the mouse data w/ human genes
obj_mmAdult = obj_mmData[,colnames(obj_adultMouse)]

# add back the metadata
obj_mmAdult = AddMetaData(obj_mmAdult, obj_adultMouse[[]])
head(obj_mmAdult[[]])
SaveH5Seurat(obj_mmAdult, here('data/tidy_data/rdas/mmData.comp.adult.raw.mm.h5seurat'))
