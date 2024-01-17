ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

suppressMessages(library(ArchR))
addArchRThreads(threads = 5)

# path to archr project
PROJDIR='/projects/pfenninggroup/singleCell/Macaque_SealDorsalHorn_snATAC-seq'
#ARCHRPROJ=paste(PROJDIR,'data/tidy_data/ArchR_Seal001', sep='/')
MJL_ARCHRPROJ = paste(PROJDIR,'data/tidy_data/ArchRProjects/Mouse_DorsalHorn_scATAC', sep='/')

# load project
proj <- loadArchRProject(path = MJL_ARCHRPROJ, force = FALSE, showLogo = FALSE)


### Notes: QC TSS vs. uniqueFrags visualized on local machine. Highest density was above 10^(3.5) frags for all 5 samples, so I selected 3.5 as the nFrags cutoff for subsequent analysis

print('cutoffs...')
print('')
print('')
## increase unique fragments cutoff to 10^3.5, remove cluster of low QC cell

idxSample <- BiocGenerics::which(proj$nFrags > 10^3.5)  # based on visual appearance of cluster in arrow QC
cellsSample <- proj$cellNames[idxSample]
proj = subsetCells(ArchRProj = proj, cellNames = cellsSample)

### previously looped over varFeats, then picked 200
varFeat = 200
# for (varFeat in c(30, 70, 110, 150, 190, 230)){
    
    print('VarFeat:')
    print(varFeat)
  # add iterative LSI
  pd = getCellColData(proj)
  dimRed = names(attributes(proj)$reducedDims)
  embedNames = names(attributes(proj)$embeddings)
  
  iterLSIName = paste0("IterativeLSI",varFeat)
  print(iterLSIName)
  if (iterLSIName %ni% dimRed){
    pdf()
  proj <- addIterativeLSI( proj, useMatrix = "TileMatrix", 
    name = iterLSIName,
    LSIMethod = 2, 
    iterations = 4, # increase this if noticing subtle batch effects
    scaleTo = 20000, # median unique fragment per cell
    selectionMethod = 'var',
    clusterParams = list( # See Seurat::FindClusters
      resolution = .2, # lower this if noticing subtle batch effects
      sampleCells = 10000,  n.start = 10), 
    varFeatures = varFeat * 1000, # also can reduce this if noticing subtle batch effects
    dimsToUse = 1:30, force = FALSE)
  dev.off()
  proj = saveArchRProject(ArchRProj = proj)}
  
  # add umap
  UMAPName = paste0("UMAPI",varFeat)
  if (UMAPName %ni% embedNames){
    print(UMAPName)
    proj <- addUMAP(proj, reducedDims = iterLSIName, 
                     name = UMAPName, nNeighbors = 30, minDist = 0.5, 
                     metric = "cosine", force = F)}

  # add clusters
  ClustersName = paste0("ClustersI",varFeat)
  if (ClustersName %ni% names(pd)){
    print(ClustersName)
    proj <- addClusters(proj, reducedDims = iterLSIName, method = "Seurat", 
                      filterBias = TRUE, name = ClustersName, resolution = .5, force = F)}
  
  # add Harmony batch correction
  HarmonyName = paste0("HarmonyI",varFeat)
  if (HarmonyName %ni% dimRed ){
    print(HarmonyName)
  proj <- addHarmony(proj, reducedDims = iterLSIName, 
                      max.iter.harmony = 15, name = HarmonyName, 
                      groupBy = c("Sample","Biological.rep"), force = F)}
  
  # add umap
  UMAPName2 = paste0("UMAPH",varFeat)
  if (UMAPName2 %ni% embedNames){
    print(UMAPName2)
  proj <- addUMAP(proj, reducedDims = HarmonyName, 
                   name = UMAPName2, nNeighbors = 30, minDist = 0.5, 
                   metric = "cosine", force = F)}
  
  # add clusters
  ClustersName2 = paste0("ClustersH",varFeat)
  if (ClustersName2 %ni% names(pd)){
    print(ClustersName2)
  proj <- addClusters(proj, reducedDims = HarmonyName, method = "Seurat", 
                      name = ClustersName2, resolution = .5, force = F)}
  proj = saveArchRProject(ArchRProj = proj)
# }

print('Done - no error')
