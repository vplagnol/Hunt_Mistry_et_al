

##default parameters
choice.phenos <- c('control', 'allAID')
code <- 'data/genotypes/NCF2_full'
#code <- 'data/genotypes/PTPRK_full'

source('scripts/association/intro_parse_data.R')


#### define output files
output.folder.calpha <- 'data/pvalues_all_samples_calpha'
output.file.calpha <- paste(output.folder.calpha, '/', base.code, '_', pheno.code, '.out', sep = '')
output.table.calpha <- paste(output.folder.calpha, '/', base.code, '_', pheno.code, '.tab', sep = '')

if (!file.exists(output.folder.calpha)) dir.create(output.folder.calpha)

output.file.single <- paste('data/pvalues_all_samples_single/', base.code, '_', pheno.code, '.out', sep = '')


output.file.burden <- paste('data/pvalues_all_samples_burden/', base.code, '_', pheno.code, '.out', sep = '')
output.file.burden.regression <- paste('data/pvalues_all_samples_burden_regression/', base.code, '_', pheno.code, '.out', sep = '')




######### define the main data frame
genotypes.full <- t(genotypes)  ##I prefer like this
pheno.table <- read.table(pheno.files[1])
pheno.values <- pheno.table$V2[ match(dimnames(genotypes.full)[[1]], table = pheno.table$V1) ] -1
my.frame <- data.frame(pheno = pheno.values)

good <- !is.na(my.frame$pheno)
genotypes.full <- genotypes.full[ good,]
my.frame <- my.frame[ good, , drop = FALSE]


### Now after removing bad samples, one can remove monomorphic SNPs
poly.num <- as.numeric(apply(genotypes.full, MAR = 2, FUN = sum, na.rm = TRUE))
poly <- poly.num > 0
genotypes.full <- genotypes.full[, poly]
summ <- summ[ poly, ]

############ Now the unconditional burden test
genotypes.cand <- genotypes.full[, summ$candidate ]
sum.cand.variants <-   apply(genotypes.cand, FUN = sum, MAR = 1, na.rm = TRUE)
MAF.candidates <- signif(sum(sum.cand.variants)/ (2*nrow(my.frame)), 3)
sum.cand.variants.bin <- sum.cand.variants > 0

my.test <- fisher.test(table(sum.cand.variants.bin, my.frame$pheno))
my.p.burden <- signif(my.test$p.value, 3)
OR <- signif(my.test$estimate, 3)
OR.CI <- paste(signif(my.test$conf.int, 3), collapse = '-')

cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', MAF.candidates, '\t',  OR, '\t', OR.CI, '\t', my.p.burden, '\n', file = output.file.burden)
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', MAF.candidates, '\t',  OR, '\t', OR.CI, '\t', my.p.burden, '\n')

mod <- glm (my.frame$pheno ~ sum.cand.variants, family = binomial)

my.coeffs <- coef(summary(mod))
OR <- signif(exp(my.coeffs[2,1]), 3)
OR.CI <- signif( exp(c(my.coeffs[2,1] - 1.96*my.coeffs[2,2], my.coeffs[2,1] + 1.96*my.coeffs[2,2])), 3)
OR.CI <- paste(OR.CI, collapse = '-')

logOR.sd <- signif(my.coeffs[2,2], 3)

my.p.burden <- signif(my.coeffs[2,4], 3)

cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', MAF.candidates, '\t',  OR, '\t', logOR.sd, '\t', OR.CI, '\t', my.p.burden, '\n', file = output.file.burden.regression)
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', MAF.candidates, '\t',  OR, '\t', logOR.sd, '\t', OR.CI, '\t', my.p.burden, '\n')
#stop()

############################################################################# Now the unconditional single SNP analysis
if (file.exists(output.file.single )) file.remove(output.file.single)
my.formula <- 'pheno ~ count.rare.variants'
for (col in 1:ncol(genotypes.full)) {
  my.count <- sum( as.numeric(genotypes.full[, col]), na.rm = TRUE)

  
  if (my.count > 0) {
    my.frame$count.rare.variants <- as.numeric(genotypes.full[, col])



    
    my.alleles <- subset(my.frame$count.rare.variants, my.frame$pheno == 0)
    my.maf <- mean(my.alleles, na.rm = TRUE)/2
    geno.controls <- pmin(subset(my.frame, pheno == 0 & !is.na(count.rare.variants), 'count.rare.variants', drop = TRUE), 1)
    geno.cases <- pmin(subset(my.frame, pheno == 1 & !is.na(count.rare.variants), 'count.rare.variants', drop = TRUE), 1)
    
    my.mat <- matrix(data = c(2*length(geno.controls) - sum(geno.controls), sum(geno.controls), 2*length(geno.cases) - sum(geno.cases), sum(geno.cases)),
                     nrow = 2,
                     ncol = 2)
    
    my.test <- fisher.test(my.mat)
    OR.CI <- paste(signif(my.test$conf.int, 3), collapse = '-')
    OR <- signif(my.test$estimate, 3)    
    my.p.SNP <- signif(my.test$p.value, 3)

                                        #my.mod <- glm (my.frame$pheno ~ my.frame$count.rare.variants, family = binomial)
                                        #my.coeffs <- coef(summary(my.mod))
                                        #my.p.GLM <- signif(my.coeffs[2,4], 3)
    
    my.name <- dimnames(genotypes.full)[[2]][ col ]
    ##cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', my.maf, '\t', my.count, '\t', my.name, '\t', OR, '\t', OR.CI, '\t', my.p.SNP, '\t', my.P.GLM, '\n', file = output.file.single, append = TRUE)
    ##cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', my.maf, '\t', my.count, '\t', my.name, '\t', OR, '\t', OR.CI, '\t', my.p.SNP, '\t', my.P.GLM, '\n')

    cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', my.maf, '\t', my.count, '\t', my.name, '\t', OR, '\t', OR.CI, '\t', my.p.SNP,  '\n', file = output.file.single, append = TRUE)
    cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', my.maf, '\t', my.count, '\t', my.name, '\t', OR, '\t', OR.CI, '\t', my.p.SNP, '\n')
  }
}

message("Output in ", output.file.single)
stop()

############################################################################# And now the C-alpha test
message('Now the c-alpha test')
#if (file.exists(output.file.calpha )) file.remove(output.file.calpha)
#if (file.exists(output.table.calpha )) file.remove(output.table.calpha)


diallelic <- summ$nalleles == 1 & summ$candidate
geno.di <- t(genotypes.full[ , diallelic ])
summ.di <- summ[ diallelic, ]

### compute the c-alpha statistic
pheno.full.loc <- factor( choice.phenos[ my.frame$pheno + 1 ], levels = choice.phenos)
proportions <-  table(pheno.full.loc)/length(pheno.full.loc)

source('scripts/association/calpha_routines.R')
my.counts <- get.count.table (geno.di, pheno = pheno.full.loc, summ.di$signature)
my.T <- calpha.test( my.counts, proportions = proportions, phenos = choice.phenos, verbose = FALSE)
message('True value is ', signif(my.T, 3))

cat(proportions, '\n', sep = '\t',  file = output.table.calpha)
write.table (x = my.counts, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = '\t', append = TRUE, file = output.table.calpha)
print(output.table.calpha)


##### Now simulate extensively
nsim <- 10000
my.T.sim <- c()

for (j in 1:nsim) {
  set.seed(j)
  my.counts <- get.count.table (geno.di, pheno = sample(pheno.full.loc, replace = FALSE), summ.di$signature)
  my.T.sim <- c(my.T.sim, calpha.test( my.counts, proportions = proportions[ choice.phenos ], phenos = choice.phenos))
  message(j, ' ', my.T.sim[j])
}


#### print it all
my.V <- var(my.T.sim)
my.chisq <- my.T^2/ my.V
my.p <- pchisq(q = my.chisq, df = 1, lower.tail = FALSE)
my.p.safe <- sum( my.T.sim > my.T)/length(my.T.sim)

#cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', signif(my.p.safe, 3), '\t', signif(my.p, 3), '\n', file = output.file.calpha)
#cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', signif(my.p.safe, 3), '\t', signif(my.p, 3), '\n')

cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', signif(my.p.safe, 3), '\n', file = output.file.calpha)
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', signif(my.p.safe, 3), '\n')


