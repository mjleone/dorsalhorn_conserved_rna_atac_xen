### conda activate r4

library(tidyverse)
library(data.table)
library(RColorBrewer)
library(here)
library(readr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#DATADIR='data/raw_data/gwas_enrichments'
main_DATADIR='data/tidy_data/gwas_enrichments_mouse'
main_CODEDIR='code/raw_code/gwas_enrichments_mouse'
main_FIGDIR='figures/exploratory/gwas_enrichments_mouse'

genebased_DATADIR='data/tidy_data/genebased_gwas_enrichments_mouse'
genebased_CODEDIR='code/raw_code/genebased_gwas_enrichments_mouse'
genebased_FIGDIR='figures/exploratory/genebased_gwas_enrichments_mouse'

liver_DATADIR='data/tidy_data/gwas_enrichments_mouse_pfenliver'
liver_CODEDIR='code/raw_code/gwas_enrichments_mouse_pfenliver'
liver_FIGDIR='figures/exploratory/gwas_enrichments_mouse_pfenliver'

immune_DATADIR='data/tidy_data/gwas_enrichments_immunecalvin'
immune_CODEDIR='code/raw_code/gwas_enrichments_immunecalvin'
immune_FIGDIR='figures/exploratory/gwas_enrichments_immunecalvin'

catlas_DATADIR='data/tidy_data/gwas_enrichments_mouse_catlas'
catlas_CODEDIR='code/raw_code/gwas_enrichments_mouse_catlas'
catlas_FIGDIR='figures/exploratory/gwas_enrichments_mouse_catlas'

## pain list
pain_list = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/deepika_pain_gwas/pain_gwas_list.txt'
pain_df = read.csv(pain_list, header = F)
pain_files = pain_df[,1]
pain_basenames = basename(pain_files)
pain_list_traits = sapply(strsplit(pain_basenames,"-"), '[', 1)



### Get Johnston stats
### don't re-run -- uses a lot of RAM
if (FALSE){
JOHNSTONDIR = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Johnston_chronic_pain-bgen.stats.gz'
JOHNSTON_SNP = read_tsv(JOHNSTONDIR)
JOHNSTON_sig_snps = sum(JOHNSTON_SNP$P_BOLT_LMM < 0.05)
print(JOHNSTON_sig_snps)
print('done')
}

### Get Kupari stats
### don't re-run -- uses a lot of RAM
if (FALSE){
KUPARIDIR = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Kupari_et_al_chronic_pain_GWAS_UKBB'
KUPARI_SNP_files = list.files(KUPARIDIR, pattern = "stats")
KUPARI_SNP_files
KUPARI_SNP_list = lapply(file.path(KUPARIDIR, KUPARI_SNP_files), read_tsv)
KUPARI_sig_snps = lapply(KUPARI_SNP_list, function(df){sum(df$P_BOLT_LMM < 0.05)})
}

# [1] "stats_all_QC.f2956_chronic_MF.txt.gz" --- General pain for 3+ months
# [2] "stats_all_QC.f3404_chronic_MF.txt.gz" --- Neck/shoulder pain for 3+ months
# [3] "stats_all_QC.f3414_chronic_MF.txt.gz" --- Hip pain for 3+ months
# [4] "stats_all_QC.f3571_chronic_MF.txt.gz" --- Back pain for 3+ months
# [5] "stats_all_QC.f3741_chronic_MF.txt.gz" --- Stomach/abdominal pain for 3+ months
# [6] "stats_all_QC.f3773_chronic_MF.txt.gz" --- Knee pain for 3+ months
# [7] "stats_all_QC.f3799_chronic_MF.txt.gz" --- Headaches for 3+ months
# [8] "stats_all_QC.f4067_chronic_MF.txt.gz" --- Facial pains for 3+ months
# [9] "stats_all_QC.fSITE_chronic_MF.txt.gz" --- Number of chronic pain sites for 3+ months

# [[1]]
# [1] 504872
# [[2]]
# [1] 693612
# [[3]]
# [1] 630020
# [[4]]
# [1] 734829
# [[5]]
# [1] 570918
# [[6]]
# [1] 714486
# [[7]]
# [1] 654940
# [[8]]
# [1] 449692
# [[9]]
# [1] 783122

### get al khoury stats
### don't re-run -- uses a lot of RAM
if (FALSE){
KHOURYDIR = '/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Khoury_et_al_chronicOverlapingPain_GWAS_UKBB'
KHOURY_SNP_files = list.files(KHOURYDIR, pattern = "stats")
KHOURY_SNP_files
KHOURY_SNP_list = lapply(file.path(KHOURYDIR, KHOURY_SNP_files), read_tsv)
KOURY_sig_snps = lapply(KHOURY_SNP_list, function(df){sum(df$P_BOLT_LMM < 0.05)})
}
# [1] "stats_all_QC.1vs0_chronic_MF.txt.gz" 
# [2] "stats_all_QC.2+vs0_chronic_MF.txt.gz"
# [3] "stats_all_QC.2+vs1_chronic_MF.txt.gz"

# [1] 838179
# [[2]]
# [1] 1112818
# [[3]]
# [1] 503343


# traits for combined figure
ordered_traits_clearer = c("Lean Body Mass", "Bone Mineral Density", "Coronary Artery Disease",
"General Pain (J)", "Local Pain (Kh)", "Multisite Pain (Kh)", "Multisite vs Local Pain (Kh)",
"Number of Pain Sites (Kh)", "General Pain (Ku)", "Headache (Ku)","Facial Pain (Ku)", "Neck/Shoulder Pain (Ku)", "Back Pain (Ku)", "Hip Pain (Ku)", "Stomach/Abdominal Pain (Ku)", "Knee Pain (Ku)")

ordered_sig_snps = c("","","",
  "","838179","1112818", "503343", 
  "783122", "504872", "654940", "449692", "693612", "734829", "630020", "570918", "714486" )

trait_clearer = c("Bone Mineral Density", "Coronary Artery Disease","Back Pain (Ku)",
                   "General Pain (J)", "Local Pain (Kh)", "Multisite Pain (Kh)", "Multisite vs Local Pain (Kh)", "Facial Pain (Ku)", "General Pain (Ku)", "Headache (Ku)", "Hip Pain (Ku)", 
                   "Knee Pain (Ku)", "Lean Body Mass", "Neck/Shoulder Pain (Ku)", "Number of Pain Sites (Kh)", "Stomach/Abdominal Pain (Ku)")

#trait_labels = paste(trait_clearer, '-', ordered_sig_snps[match(trait_clearer, ordered_traits_clearer)])
#trait_labels= paste(ordered_traits_clearer, '-', ordered_sig_snps)
trait_labels = ordered_traits_clearer

# traits_coded = c("BMD", "CAD","BackPain", 
#                    "ChronicPain", "CLPC", "COPC", "COPCvCLPC","FacialPain", "GeneralPain", "Headache","HipPain",
#                    "KneePain", "LBM", "NeckShoulderPain", "NumCPSites", "StomachAbdominalPain")


# # just traits in the figure                   
# pain_files = pain_files[pain_list_traits %in% traits_coded]
# pain_list_traits = pain_list_traits[pain_list_traits %in% traits_coded]

# # get the right order or files, and extract sample sizes
# pain_files_ordered =pain_files[match(traits_coded,pain_list_traits) ]
# sample_sizes = lapply(pain_files_ordered, function(filename){

#   df = read.table(gzfile(filename), header = TRUE, fill = TRUE)

#   return(sum(abs(df$Z[!is.na(df$Z)]) > 1.645)  )   })
  # # if there's just one sample size, return
  # if ( length(unique(df$N[!is.na(df$N)])) == 1){
  #   N = unique(df$N[!is.na(df$N)])[1]
  #   return(N)
  # }
  # else{
  #   return(NA)
  # }

  # })


# celltypes of spinal cord

new_neuron_names = c('Exc-BNC2/HMGA2','Exc-NMUR2/PREX2', 'Exc-NMU/TAC3', 'Exc-TAC3/COL5A2/PLCH1','Exc-LMO3/TRPC3', 'Exc-SKOR2/NELL2', 'Exc-MAF/ADARB2', 'Exc-MAFA/BNC2', 
                  'Exc-NTS/TSHZ2', 'Exc-SNTB1/TRH/DACH1', 'Exc-PBX3/SLIT2', 
                  'Inh-CACNA2D3/TCF4', 'Inh-SORCS1/SDK1', 'Inh-PDZD2/SGCD', 'Inh-NPY/ZIC1', 'Inh-MEF2C/RORB','Inh-NXPH1/SDK2','Inh-PDYN/PTPRK')  

reverse_celltype_order = c("Oligo 1",
                 "Oligo 2", 
                 "Astrocyte 1",
                 "Astrocyte 2",
                 "OPC",
                 "Endothelial",
                 "Ependymal",
                 "Microglia",
                 "Mural",
                 "Meninges",
                 " ",
                 new_neuron_names)
                 
celltypes = c(paste0('GLUT', 1:11),
              'GABA1', 'GABA2_1','GABA2_2','GABA3',
            'GABA4_1','GABA4_2', 'GABA5',
                 "Oligo.1",
                 "Oligo.2", 
                 "Astrocyte.1",
                 "Astrocyte.2",
                 "OPC",
                 "Endothelial",
                 "Ependymal.cells",
                 "Microglia",
                 "Mural",
                 "Meninges")
               
new_celltypes = c(new_neuron_names,
                 "Oligo 1",
                 "Oligo 2", 
                 "Astrocyte 1",
                 "Astrocyte 2",
                 "OPC",
                 "Endothelial",
                 "Ependymal",
                 "Microglia",
                 "Mural",
                 "Meninges")
                 
names(new_celltypes) = celltypes

# color accordingly 
celltypes_col = c(brewer.pal(11, 'Paired'), 
                  brewer.pal(7, 'Dark2'), 
                  brewer.pal(10, 'Set3'))


# add names vector in order to name the colors
names(celltypes_col) = celltypes


### "UCLA_NAIVE.idr.optimal_peak"                 
immune_cells = c("fullard_putaman_optimal_peak", "fullard_hippocampus_neurons.idr.optimal_peak", 
"UCLA_NAIVE.idr.optimal_peak","UCLA_IFNB.idr.optimal_peak", "UCLA_IFNG.idr.optimal_peak")
immune_cells_clearer = c("putaman neurons", "hippocampal neurons", 
"Naive immune cells", "IFNB macrophages", "IFNG macrophages")

names(immune_cells_clearer) = immune_cells


###
###
###
##

enrich_target_neg = readRDS("enrich_target_neg.rds")
genebased_enrich_target_neg = readRDS("genebased_enrich_target_neg.rds")
genebased_enrich_target_neg = genebased_enrich_target_neg %>% filter(class != 'GLIA')
immune_enrich_target_neg = readRDS("immune_enrich_target_neg.rds")
immune_enrich_target_neg$new_celltypes = immune_cells_clearer[immune_enrich_target_neg$celltype]

liver_enrich_target_neg = readRDS("liver_enrich_target_neg.rds")
liver_enrich_target_neg$new_celltypes = "bulk liver"

pdf('just_spinal_cord_test.pdf', width = 7,  height = 4)

ggplot(enrich_target_neg, aes(x = factor(new_celltypes, levels= reverse_celltype_order), 
y = factor(trait_clearer, levels= ordered_traits_clearer), fill = -log10(FDR))) +
  geom_tile(data=enrich_target_neg, aes(color=signif.FDR, width=.9, height=.9), size=1) +
  scale_fill_viridis_c() +
  scale_fill_gradient2(limits=c(0, 6)) +
  scale_color_manual("FDR < 0.05", values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1), text = element_text(size=10)) +
  scale_y_discrete(limits=rev) +
  labs(x = "Cell Type", y = "GWAS Trait")

dev.off()


## remove GLUT9, 11, and GABA2_1
to_remove = c('Exc-NTS/TSHZ2', 'Exc-PBX3/SLIT2', 'Inh-SORCS1/SDK1')
new_neuron_names_removed = new_neuron_names[!(new_neuron_names %in% to_remove)]
genebased_enrich_target_neg = genebased_enrich_target_neg[!(genebased_enrich_target_neg$new_celltypes %in% to_remove),]


pdf('just_genebased_test.pdf', width = 4.41,  height = 3)

ggplot(genebased_enrich_target_neg, aes(x = factor(new_celltypes, levels= new_neuron_names_removed), 
y = factor(trait_clearer, levels= ordered_traits_clearer), fill = -log10(FDR))) +
  geom_tile(data=genebased_enrich_target_neg, size=.05, color = 'grey70') +
  scale_fill_viridis_c() +
  scale_fill_gradient2(limits=c(0, 6)) +
  scale_color_manual("FDR < 0.05", values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = "none", text = element_text(size=10),
  panel.border = element_rect(colour = "black", fill=NA, size=.5)) +
  scale_y_discrete(limits=rev, labels = rev(trait_labels) )  +
  labs(x = "", y = "")

dev.off()


#### combining
full_cell_order = c("bulk liver", immune_cells_clearer, reverse_celltype_order)
enrich_combined = rbind(liver_enrich_target_neg, enrich_target_neg, immune_enrich_target_neg, fill = TRUE)

## remove GLUT9, 11, and GABA2_1
to_remove = c('Exc-NTS/TSHZ2', 'Exc-PBX3/SLIT2', 'Inh-SORCS1/SDK1')
new_neuron_names_removed = new_neuron_names[!(new_neuron_names %in% to_remove)]
enrich_combined = enrich_combined[!(enrich_combined$new_celltypes %in% to_remove),]

neuron_nicknames_removed = c('Exc-BNC2','Exc-NMUR2', 'Exc-NMU', 'Exc-TAC3','Exc-LMO3', 'Exc-SKOR2', 'Exc-MAF', 'Exc-MAFA', 
                  'Exc-SNTB1', 
                  'Inh-CACNA2D3','Inh-PDZD2', 'Inh-NPY', 'Inh-MEF2C','Inh-NXPH1','Inh-PDYN')  

reverse_celltype_nickname = c("Oligo 1",
                 "Oligo 2", 
                 "Astrocyte 1",
                 "Astrocyte 2",
                 "OPC",
                 "Endothelial",
                 "Ependymal",
                 "Microglia",
                 "Mural",
                 "Meninges",
                 " ", neuron_nicknames_removed)
full_nickname_order = c("bulk liver", immune_cells_clearer,  reverse_celltype_nickname)

enrich_combined$new_nicknames = sapply(strsplit(enrich_combined$new_celltypes,"/"), '[', 1)



print(table(enrich_combined$signif.FDR))

pdf('combined_test.pdf', width = 8,  height = 4)

ggplot(enrich_combined, aes(x = factor(new_nicknames, levels= full_nickname_order), 
y = factor(trait_clearer, levels= ordered_traits_clearer), fill = -log10(FDR))) +
  geom_tile(data=enrich_combined, size=.05, color = 'grey70') +
  geom_tile(data=enrich_combined[enrich_combined$signif.FDR == TRUE,], 
    aes(color=signif.FDR, width=.95, height=.95), size=.5) + 
  scale_fill_gradient2(limits=c(0, 6)) +
  scale_color_manual("FDR < 0.05", values = c('TRUE' = 'black', 'FALSE' = 'white')) +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1), text = element_text(size=10),
  panel.border = element_rect(colour = "black", fill=NA, size=.5), legend.box.spacing = unit(6, "pt"), legend.margin=margin(0,0,0,0)) +
  scale_y_discrete(limits=rev, labels = rev(trait_labels) )  +
  labs(x = "Cell Type", y = "GWAS Trait")

  # name cell types
  #scale_x_discrete(labels= nicknames)

dev.off()
