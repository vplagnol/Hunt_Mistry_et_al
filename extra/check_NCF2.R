load('data/genotypes/NCF2_full.RData')

coeliac <- read.table('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/phenotypes/coeliac.phe')
allAID <- read.table('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/phenotypes/allAID.phe')


allAID.v2 <- subset(allAID, ! V1 %in% subset(coeliac, V2 == 2, 'V1', drop = TRUE))


SNP1 <- which(summ$dbSNP137 == 'rs17849502')
SNP2 <- which(summ$dbSNP137 == 'rs17849501')


geno1 <- genotypes[SNP2, dimnames(genotypes)[[2]] %in% allAID.v2$V1 ]
pheno1 <- allAID$V2[ match(table = allAID$V1, names(geno1)) ]
tab1 <- table(geno1, pheno1)
fisher.test(rbind(tab1[1,]*2 + tab1[2,], tab1[2,] + 2*tab1[3,]))
