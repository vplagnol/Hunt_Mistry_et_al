library(snpStats)

region <- 'data/genotypes/IL2_full'

load(paste(region, '.RData', sep = ''))

ichip <- read.plink('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/fluidigm50k_final15jan2013_noIFCbarcode.i-view_indiv_withphenotype_plink__toload') 