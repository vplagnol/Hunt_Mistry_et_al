source('scripts/association/qq_chisq_VP.R')
#library(snpStats)


##############################
message('Unconditional burden analysis')
system("cat /data_n2/vplagnol/Projects/fluidigm/association_v2/data/pvalues_all_samples_burden/*out > temp/all_burden.tab")

burden.uncond <- read.table('temp/all_burden.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'MAF.candidates', 'OR.CI', 'OR', 'P.burden'))
write.csv(x = burden.uncond, file = 'results/combined_tables_supplemental/burden_unconditional_fishertest.csv', row.names = FALSE)


##############################
message('Conditional burden analysis with ichip')
system("cat /data_n2/vplagnol/Projects/fluidigm/association_v2/data/pvalues_with_ichip/*out > temp/all_P_ichip.tab")

data <- read.table('temp/all_P_ichip.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'count', 'ichip.covariates', 'ichip.P', 'rare.MAF', 'rare.P'))
write.csv(x = data, file = 'results/combined_tables_supplemental/burden_conditional_glm.csv', row.names = FALSE)


#my.res.matrix <- matrix( rep(NA, nrow(data)), dimnames = list(unique(data$locus), unique(data$phenos)), nrow = 25, ncol = 7)
my.res.matrix <- matrix( rep(NA, 25*6), dimnames = list(unique(data$locus), unique(data$phenos)), nrow = 25, ncol = 6)

for (i in 1:nrow(my.res.matrix)) {
  my.locus <- dimnames(my.res.matrix)[[1]][i]
  for (j in 1:ncol(my.res.matrix)) {
    my.phenos <- dimnames(my.res.matrix)[[2]][j]
    if (sum(data$locus == my.locus &  data$phenos == my.phenos)  == 1) {
      my.res.matrix[i,j] <- subset(data, locus == my.locus & phenos == my.phenos, 'rare.P', drop = TRUE)
    } else {
      message('Issue with ', my.phenos, ' ', my.locus)
    }
  }
}

output.tab <- 'results/pvalues_with_ichip_covariates.csv'
write.csv(x = my.res.matrix, file = output.tab)
print(output.tab)



###################
message('UNIQ analysis')
system("cat /data_n2/vplagnol/Projects/fluidigm/association_v2/data/pvalues_all_samples_UNIQ/*out > temp/all_P_UNIQ.tab")
uniq <- read.table('temp/all_P_UNIQ.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'nsim', 'uniq.p.cases', 'uniq.p.controls'))
write.csv(x = uniq, file = 'results/combined_tables_supplemental/UNIQ_permutation_test.csv', row.names = FALSE)


###################
message('Single SNP analysis with maximum sample size')
system('cat data/pvalues_all_samples_single/* | sort -gk10 > temp/pvalues_single_all_samples.tab')
data2 <- read.table('temp/pvalues_single_all_samples.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'rare.maf', 'count', 'rare.name', 'OR', 'OR.CI', 'rare.P'))
data2 <- data2[ order(data2$rare.P), ]

output.tab <- 'results/pvalues_single_all_samples.tab'
write.csv(x = data2, file = output.tab, row.names = FALSE)
print(output.tab)
write.csv(x = data2, file = 'results/combined_tables_supplemental/pvalues_single_sample_unconditional_fishertest.csv', row.names = FALSE)


###########
system('cat data/pvalues_all_samples_calpha/*out | sort -rgk5 >  temp/pvalues_calpha_all_samples.tab')
data3 <- read.table('temp/pvalues_calpha_all_samples.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'p.value'))
data3 <- data3[ order(data3$p.value), ]
write.csv(x = data3, file = 'results/combined_tables_supplemental/pvalues_calpha.csv', row.names = FALSE)


output.tab <- 'results/pvalues_calpha_all_samples.tab'
write.csv(x = data3, file = output.tab, row.names = FALSE)
print(output.tab)







#############################################################################

system('cat data/pvalues_with_ichip_single/* | sort -rgk9 > temp/pvalues_single_with_ichip.tab')
data1 <- read.table('temp/pvalues_single_with_ichip.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'count', 'ichip.covariates', 'ichip.P', 'rare.maf', 'rare.name', 'rare.P'))
data1 <- data1[ order(data1$rare.P), ]
write.csv(x = data1, file = 'results/combined_tables_supplemental/pvalues_singleSNP_conditionaliChip_glm.csv', row.names = FALSE)



#output.pdf <- 'qqplot_t1d_coeliac_single_at_least_four_calls.pdf'
#pdf(output.pdf,width = 8, height = 4)
#par(mfrow = c(1, 2))
#p.coeliac <-  subset(data1, phenos == 'control_coeliac', 'rare.P', drop = TRUE) 
#qq.chisq(x = -2*log(p.coeliac), df = 2, x.max = 30, main = 'Coeliac: coding variants p-values', pch = '+', slope.lambda = TRUE, overdisp = TRUE)
#p.t1d <-  subset(data1, phenos == 'control_t1d',  'rare.P', drop = TRUE) 
#qq.chisq(x = -2*log(p.t1d), df = 2, x.max = 30, ylim = c(0, 30), main = 'T1D: coding variants p-values', pch = '+', slope.lambda = TRUE, overdisp = TRUE)
#dev.off()
#print(output.pdf)




output.tab <- 'results/pvalues_single_with_ichip_covariates.tab'
write.csv(x = data1, file = output.tab, row.names = FALSE)
print(output.tab)



################ create Figure 1 with allAID data



#my.conf <-  c(0.025, 0.975)
my.conf <-  c(0.005, 0.995)

output.pdf <- 'figure1.pdf'
pdf(output.pdf,width = 3, height = 13)
par(mfrow = c(5, 1))




my.trim <- 0.6
data2 <- subset(data2, rare.maf < 0.051 & count > 6)  ##unconditional
data1 <- subset(data1, rare.maf < 0.051 & count > 6)  ##conditional

#p.allAID <-  subset(data2, phenos == 'control_allAID', 'rare.P', drop = TRUE)  ##unconditional
p.allAID <-  subset(data1, phenos == 'control_allAID', 'rare.P', drop = TRUE)  ##conditional
message('Number of tests: ', length(p.allAID))


#for (i in 1:2) {
  #if (i == 1) pdf('fig/figure1a.pdf')
  #qq.chisq.VP(x = -2*log(p.allAID), ##1- conditional single-SNP P-values
  #            df = 2,
  #            x.max = 30,
  #            trim = my.trim,
  #            main = 'a- All AID, single variant tests',
  #            pch = '+',
  #            slope.lambda = TRUE,
  #            ylim = c(0, 30),
  #            conc = my.conf,
  #            overdisp = FALSE,
  #            sub = paste(length(p.allAID), 'variants (MAF < 5%, > 6 non ref alleles)')) 
  #if (i == 1) dev.off()
#}

## 2- C-alpha pooled
data3 <- subset(data3, ncases > 420)

for (i in 1:2) {
  if (i == 1) pdf('fig/figure1b.pdf')
  qq.chisq.VP(x = -2*log(data3$p.value),
              df = 2,
              x.max = 30,
              ylim = c(0, 30),
              main = 'a- C-alpha test (gene based)',
              pch = '+',
              trim = my.trim,
              conc = my.conf,
              slope.lambda = TRUE,
              sub = '175 tests (25 genes by 7 cohorts)\n(functional variants with MAF < 0.5%)',
              line.sub = 4)
  if (i == 1) dev.off()
}

##3- burden unconditional
my.p <- subset(burden.uncond, MAF.candidates > 0.005, 'P.burden', drop = TRUE)
for (i in 1:2) {
  
  if (i == 1) pdf('fig/figure1c.pdf')
  qq.chisq.VP(x = -2*log(my.p),
              df = 2,
              x.max = 30,
              ylim = c(0, 30),
              main = 'b- Burden test (gene based)',
              pch = '+',
              conc = my.conf,
              slope.lambda = TRUE,
              trim = my.trim,
              sub = '175 tests (25 genes by 7 cohorts)\n(functional variants with MAF < 0.5%)',
              line.sub = 4)
  if (i == 1) dev.off()
}

## 4- burden conditional
my.p <- subset(data, count > 40, 'rare.P', drop = TRUE)
for (i in 1:2) {
  if (i == 1) pdf('fig/figure1d.pdf')
  qq.chisq.VP(x = -2*log(my.p),
              df = 2,
              x.max = 30,
              ylim = c(0, 30),
              conc = my.conf,
              main = 'c- Conditional burden test (gene based)',
              pch = '+',
              slope.lambda = TRUE,
              trim = my.trim,
              sub = '150 tests (25 genes by 6 cohorts)\n(functional variants with MAF < 0.5%)',
              line.sub = 4)
  if (i == 1) dev.off()
}


##5- UNIQ test- cases
for (i in 1:2) {
  if (i == 1) pdf('fig/figure1e.pdf')
  qq.chisq.VP(x = -2*log(uniq$uniq.p.cases),
              df = 2,
              x.max = 30,
              ylim = c(0, 30),
              main = 'd- UNIQ in cases test (gene based)',
              pch = '+',
              conc = my.conf,
              trim = my.trim,
              slope.lambda = TRUE,
              sub = '175 tests (25 genes by 7 cohorts)\n(functional variants with MAF < 0.5%)',
              line.sub = 4)
  if (i == 1) dev.off()
}

##6- UNIQ test- controls
qq.chisq.VP(x = -2*log(uniq$uniq.p.controls),
            df = 2,
            x.max = 30,
            ylim = c(0, 30),
            conc = my.conf,
            main = 'e- UNIQ in controls test (gene based)',
            pch = '+',
            trim = my.trim,
            slope.lambda = TRUE,
            sub = '175 tests (25 genes by 7 cohorts)\n(functional variants with MAF < 0.5%)',
            line.sub = 4)


dev.off()
print(output.pdf)


######## create supp Figure with all cohorts
output.pdf <- 'figureS1.pdf'
pdf(output.pdf, width = 6, height = 9)
par(mfrow = c(3, 2))


my.diseases <-  c('t1d', 'coeliac', 'AITD', 'ms', 'crohns', 'allAID')
my.titles <- c('T1D', 'Coeliac disease', 'AITD', 'Multiple sclerosis', 'Crohn disease', 'All AID')
for (i in 1:6) {

  disease <- my.diseases[i]
  my.p <-  subset(data1, phenos == paste('control_', disease, sep = ''), 'rare.P', drop = TRUE) 
  qq.chisq.VP (x = -2*log(my.p),
               df = 2,
               x.max = 30,
               conc = my.conf,
               ylim = c(0, 30),
               main = my.titles[i],
               pch = '+',
               slope.lambda = TRUE,
               overdisp = FALSE,
               trim = 0.9,
               line = 4,
               sub = 'Single SNP tests (MAF < 5%, non ref count > 6)\nConditional on iChip covariates') ##1- conditional single-SNP P-values
}
dev.off()
print(output.pdf)


coeliac.interesting <- subset(data1, phenos == 'control_coeliac' & rare.P < 0.01)
crohn.interesting <- subset(data1, phenos == 'control_crohns' & rare.P < 0.01)


coeliac.interesting <- subset(data1, phenos == 'control_coeliac' & rare.P < 0.01)
crohn.interesting <- subset(data1, phenos == 'control_crohns' & rare.P < 0.01)


write.csv(x = coeliac.interesting, row.names = FALSE, 'coeliac_interesting_singleSNP_conditional.csv')
write.csv(x = crohn.interesting, row.names = FALSE, 'crohn_interesting_singleSNP_conditional.csv')

