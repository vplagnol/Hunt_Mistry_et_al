system('cat data/pvalues_all_samples_burden_regression/* > temp/burden_regression.tab')
burden.uncond <- read.table('temp/burden_regression.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'MAF.candidates', 'OR', 'logOR.sd', 'OR.CI', 'P.burden'))
burden.uncond <- subset(burden.uncond, phenos != 'control_allAID')

burden.uncond$max.OR <- log(as.numeric(gsub(pattern = '.*-', replacement = '', burden.uncond$OR.CI)))
burden.uncond$min.OR <- log(as.numeric(gsub(pattern = '-.*', replacement = '', burden.uncond$OR.CI)))

burden.uncond <- subset(burden.uncond, logOR.sd < 10)

burden.uncond$min.OR <- pmax(log(1/50), burden.uncond$min.OR )
burden.uncond$max.OR <- pmin(log(50), burden.uncond$max.OR )

f.common <- 0.2
OR.common <- 1.2
var.common <- log(OR.common)^2*f.common*(1-f.common)


burden.uncond$var.risk <- burden.uncond$max.OR^2*burden.uncond$MAF.candidates*(1-burden.uncond$MAF.candidates)
burden.uncond$var.prot <- burden.uncond$min.OR^2*burden.uncond$MAF.candidates*(1-burden.uncond$MAF.candidates)
burden.uncond$var <- pmax(log(burden.uncond$OR), 1/50)^2*burden.uncond$MAF.candidates*(1-burden.uncond$MAF.candidates)

print(mean(burden.uncond$var/var.common))

my.sim <- c()
for (nsimul in 1:1000) {
  burden.uncond$sim.OR <- exp(rnorm(log(burden.uncond$OR), sd = burden.uncond$logOR.sd))
  
  burden.uncond$var.sim <- pmax(log(burden.uncond$sim.OR), 1/50)^2*burden.uncond$MAF.candidates*(1-burden.uncond$MAF.candidates)
  my.sim <- c(my.sim, mean(burden.uncond$var.sim/var.common))
}

print(sort(my.sim)[c(50,500, 950)])

stop()


#print(pmax(burden.uncond$var.risk, burden.uncond$var.prot)/var.common)
#print(mean(pmax(burden.uncond$var.risk, burden.uncond$var.prot)/var.common))



#print(range(pmax(burden.uncond$var.risk, burden.uncond$var.prot)/var.common))



stop()






uniq <- read.table('temp/all_P_UNIQ.tab', header = FALSE, col.names = c('locus', 'phenos', 'ncontrols', 'ncases', 'nsim', 'uniq.p.cases', 'uniq.p.controls'), stringsAsFactors = FALSE)
b
interesting <- subset(uniq, uniq.p.cases < 0.01 & phenos != 'control_allAID')

interesting <- subset(interesting, phenos == 'control_coeliac')

for (i in 1:nrow(interesting)) {

  file <- paste('data/genotypes/', interesting$locus [ i ], '.RData', sep = '')
  load(file)


  disease <- gsub(pattern = '.*_', replacement = '', interesting$phenos[i])
  #for (disease in c('t1d', 'coeliac')) {
  for (disease in c('coeliac')) {

    pheno.files <- paste('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/phenotypes/', disease, '.phe', sep = '')
    pheno <- read.table(pheno.files)

### use only the polymorphic ones
    genotypes.loc <- genotypes[, dimnames(genotypes)[[2]] %in% pheno$V1 ]
    poly.num <- as.numeric(apply(genotypes, MAR = 1, FUN = sum, na.rm = TRUE))
    genotypes.loc <- genotypes.loc[ summ$nhets <= 3 & summ$candidate & poly.num >= 1, ]
    
###
    trait <- pheno$V2 [ match(dimnames(genotypes.loc)[[2]], table = pheno$V1) ] - 1
    overall.burden <- apply(MAR = 2, FUN = sum, genotypes.loc, na.rm = TRUE)
    my.tab <- table(overall.burden, trait)
    my.mod <- glm (trait ~ overall.burden, family = binomial)
    
    stop()
    #
    #count.controls <- apply(genotypes[, trait == '1'], MAR = 1, FUN = sum, na.rm = TRUE)
    #only.in.cases <- genotypes[count.controls == 0,]
    #uniq.cases <- sum(only.in.cases, na.rm = TRUE)
    #freq <- apply(MAR = 1, only.in.cases, FUN = sum, na.rm = TRUE)
    #print(table(freq))
    

  
  }
}
