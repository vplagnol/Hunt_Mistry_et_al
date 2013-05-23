
my.names <- c('Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP6500_ALL','X1000g2012apr_ALL','dbSNP137','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++','Chr','Start','End','Ref','Obs','Call', 'junk')
annovar <- read.csv('temp/annovar/fluidigm50k_final15jan2013_noIFCbarcode__has_phenotype.indiv_4377sites_41911indiv.genome_summary.csv', col.names = my.names, stringsAsFactor = FALSE, na.string = c('NA', ''))
annovar$Chr <- gsub(annovar$Chr, pattern = 'chr', replacement = '') 
annovar$signature <- paste(annovar$Chr, annovar$Start, annovar$Ref, annovar$Obs, sep = '_')



system("cat data/pvalues_all_samples_single/* | sort -grk7 > temp/maintab.tab")
data <- read.table('temp/maintab.tab', col.names = c('locus', 'pheno', 'ncontrols', 'ncases', 'rare.maf', 'count', 'signature', 'OR', 'OR.CI', 'rare.P'))

data <- subset(data, rare.maf < 0.05 & rare.P < 1.*10^(-4))
data <- merge(data, annovar[, c('Func', 'ExonicFunc', 'Ref', 'Obs', 'signature', 'dbSNP137')], by = 'signature')
data <- subset(data, !is.na(ExonicFunc))


############ now look at conditional P-valuesd
system('cat data/pvalues_with_ichip_single/* > temp/single_ichip.tab')
ichip.cond <- read.table('temp/single_ichip.tab', col.names = c('locus', 'pheno', 'ncontrols', 'ncases', 'count', 'covar', 'covar.P', 'rare.maf', 'signature', 'pval.cond'), header = FALSE, sep = '\t')
ichip.cond$signature <- gsub(pattern = ' ', replacement = '', ichip.cond$signature)
ichip.cond$pheno <- gsub(pattern = ' ', replacement = '', ichip.cond$pheno)

data <- merge(data, ichip.cond[, c('signature', 'covar', 'covar.P', 'pval.cond', 'pheno')], by = c('pheno', 'signature'))
data$Gene <- gsub(pattern = '_full', replacement = '', data$locus)
data$in.step.reg <- ifelse (data$signature == '1_183532580_G_T', 'yes', 'no')


pretty.tab <- data[, c('dbSNP137', 'ExonicFunc', 'pheno', 'Gene', 'rare.maf', 'rare.P', 'OR', 'OR.CI', 'pval.cond')]

write.csv(x = pretty.tab, row.names = FALSE, file = 'key_table.csv')
