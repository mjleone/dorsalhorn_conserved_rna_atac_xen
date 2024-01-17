#!/bin/bash

path_to_gwas="/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Kupari_et_al_chronic_pain_GWAS_UKBB/"
path_to_munged="/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/deepika_pain_gwas/"


############################################# HEADACHE f3799 ###############################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3799_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}Headache-Kupari_2021 \
	--N 200206 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# FACIAL PAIN f4067 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f4067_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}FacialPain-Kupari_2021 \
	--N 167320 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# NECK/SHOULDER PAIN f3404 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3404_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}NeckShoulderPain-Kupari_2021 \
	--N 227758 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# STOMACH/ABDOMINAL PAIN f3741 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3741_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}StomachAbdominalPain-Kupari_2021 \
	--N 182574 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# BACK PAIN f3571 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3571_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}BackPain-Kupari_2021 \
	--N 234459 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# HIP PAIN f3414 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3414_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}HipPain-Kupari_2021 \
	--N 199459 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# KNEE PAIN f3773 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f3773_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}KneePain-Kupari_2021 \
	--N 232062 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# GENERAL PAIN f2956 ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.f2956_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}GeneralPain-Kupari_2021 \
	--N 169084 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# NUM CP SITES fSITE ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.fSITE_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--out ${path_to_munged}NumCPSites-Kupari_2021 \
	--N 408247 \
	--chunksize 500000 \
	--p P_BOLT_LMM  \
	--a1 ALLELE1 \
	--a2 ALLELE0


############################################# CHRONIC PAIN JOHNSTON ################################################

python ldsc/munge_sumstats.py \
	--sumstats /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Johnston_chronic_pain-bgen.stats.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--N 387649 \
	--N-cas 169027 \
	--a1 ALLELE1 \
	--a2 ALLELE0 \
	--frq A1FREQ \
	--p P_BOLT_LMM \
	--out ${path_to_munged}ChronicPain-Johnston_2019


path_to_gwas="/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas/unmunged_sumstats/Khoury_et_al_chronicOverlapingPain_GWAS_UKBB/"

############################################# CLPC KHOURY ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.1vs0_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--N 175769 \
	--a1 ALLELE1 \
	--a2 ALLELE0 \
	--p P_BOLT_LMM \
	--out ${path_to_munged}CLPC-Khoury_2021 \
	--chunksize 500000


############################################# COPC KHOURY ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.2+vs0_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--N 164778 \
	--a1 ALLELE1 \
	--a2 ALLELE0 \
	--p P_BOLT_LMM \
	--out ${path_to_munged}COPC-Khoury_2021 \
	--chunksize 500000 

############################################# COPC vs. CLPC KHOURY ################################################

python ldsc/munge_sumstats.py \
	--sumstats ${path_to_gwas}stats_all_QC.2+vs1_chronic_MF.txt.gz \
	--merge-alleles /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist \
	--N 164778 \
	--a1 ALLELE1 \
	--a2 ALLELE0 \
	--p P_BOLT_LMM \
	--out ${path_to_munged}COPCvCLPC-Khoury_2021 \
	--chunksize 500000 


