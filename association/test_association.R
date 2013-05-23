### combined rare variant regression + single SNP analysis, all conditional on ichip covariates

##default parameters
choice.phenos <- c('control', 'ms')
code <- 'data/genotypes/ZMIZ1_full'

source('scripts/association/intro_parse_data.R')

#### define output file
output.file <- paste('data/pvalues_with_ichip/', base.code, '_', pheno.code, '.out', sep = '')
output.file.single <- paste('data/pvalues_with_ichip_single/', base.code, '_', pheno.code, '.out', sep = '')



###start by removing very rare variants in the immunochip data
ichip$genotypes <- ichip$genotypes[, col.summary(ichip$genotypes)$MAF > 0.03  & dimnames(ichip$genotypes)[[2]] != 'ichip_hg19__1__183532580']

if (pheno.code == 'control_psoriasis') stop()


########### first step that includes the immunochip
pheno.table <- read.table(pheno.files[1])
genotypes.inter.ichip <- t(genotypes[ , dimnames(ichip$genotypes)[[1]] ])

pheno.values <- pheno.table$V2[ match(dimnames(genotypes.inter.ichip)[[1]], table = pheno.table$V1) ] -1
genotypes.inter.ichip <- genotypes.inter.ichip[!is.na(pheno.values), ]
ichip.geno <- ichip$genotypes[ !is.na(pheno.values), ]
pheno.values <- pheno.values[ !is.na(pheno.values) ]


### quick regression on the ichip data
my.frame <- data.frame(pheno = pheno.values)
row.names(my.frame) <- dimnames(genotypes.inter.ichip)[[1]]

##### count nb of rare variants
my.func <- c('frameshift insertion', 'frameshift deletion', 'nonframeshift deletion',  'nonframeshift insertion', 'nonsynonymous SNV', 'stopgain SNV', 'stoploss SNV')

#threshold <- 0.01
#summ$rare <- (is.na(summ$X1000g2012apr_ALL) | summ$X1000g2012apr_ALL < threshold) &  (is.na(summ$ESP6500) | summ$ESP6500 < threshold)
#summ$candidate.v3 <- summ$rare & (  (summ$ExonicFunc %in% my.func) | (summ$Func %in% c('splicing', 'exonic;splicing')) )

genotypes.cand <- genotypes.inter.ichip[ , summ$candidate, drop = FALSE]
my.frame$count.rare.variants <- apply(genotypes.cand, MAR = 1, FUN = sum, na.rm = TRUE)
#my.frame$count.rare.variants <- rbinom(size = 2, p = 0.01, n = nrow(my.frame))

####### Do we need immunochip covariates?
covar <- c()

min.P <- 0
my.formula <- 'pheno ~ 1'
P.step <- c()
P.threshold <- 10^(-6)

while (min.P < P.threshold) {
  message(my.formula)
  my.test <- snp.rhs.tests(snp.data = ichip.geno, data = my.frame, formula = as.formula(my.formula), family = 'binomial')
  min.SNP <- which.min(p.value(my.test))
  min.P <- p.value(my.test)[ min.SNP ]
  name.SNP <- names(my.test)[ min.SNP ]
  
  if (min.P < P.threshold) {
    my.frame[, name.SNP ] <- as(as(ichip.geno[, min.SNP], 'numeric'), 'numeric')
    my.formula <- paste(my.formula, name.SNP, sep = ' + ')
    covar <- c(covar, name.SNP)
    P.step <- c(P.step, min.P)
    P.threshold <- 10^(-4)
  }
}

my.formula <- gsub(my.formula, pattern = '~ 1', replacement = '~ count.rare.variants')
#my.formula <- 'pheno ~ count.rare.variants'; covar <- c()
print(my.formula)

my.count <- sum(my.frame$count.rare.variants > 0)
#### Now do the association testing
my.test <- glm (data = my.frame, formula = my.formula, family = binomial)
my.p.rare <- coef(summary(my.test))['count.rare.variants', 4]

if (length(covar) > 0) {
  my.covar.p <- coef(summary(my.test))[covar, 4]
  covar.string <- paste(covar, collapse = '+')
  covar.p.string <- paste(signif(my.covar.p, 3), collapse = '+')
} else {
  covar.string <- 'none'
  covar.p.string <- NA
}

my.rare.MAF <- mean(my.frame$count.rare.variants, na.rm = TRUE)/2
#my.p.rare <- fisher.test(table(my.frame$count.rare.variants > 0, my.frame$pheno))$p.value

######## Now display the whole thing
cat(base.code, '\t', pheno.code, '\t', sum(pheno.values == 0), '\t', sum(pheno.values == 1), '\t', my.count,  '\t', covar.string, '\t', covar.p.string, '\t', my.rare.MAF, '\t', signif(my.p.rare, 3), '\n', file = output.file)
cat(base.code, '\t', pheno.code, '\t', sum(pheno.values == 0), '\t', sum(pheno.values == 1), '\t', my.count,  '\t', covar.string, '\t', covar.p.string, '\t', my.rare.MAF, '\t', signif(my.p.rare, 3), '\n')
message("Output in ", output.file)


############################################################################# Now the single SNP analysis, still conditional on top SNPs
if (file.exists(output.file.single )) file.remove(output.file.single) ##important to remove the old values

res.data <- data.frame(SNP = dimnames(genotypes.inter.ichip)[[2]],
                       pval.cond = NA,
                       MAF = NA)


my.counts <- apply(genotypes.inter.ichip, MAR = 2, FUN = sum, na.rm = TRUE)
genotypes.inter.ichip2 <- genotypes.inter.ichip[, my.counts > 0 ]

for (col in 1:ncol(genotypes.inter.ichip2)) {
  my.count <- sum( as.numeric(genotypes.inter.ichip2[, col]), na.rm = TRUE)
  message('Number of non ref alleles ', my.count)
  
  my.frame$count.rare.variants <- as.numeric(genotypes.inter.ichip2[, col])
  
  my.alleles <- subset(my.frame$count.rare.variants, my.frame$pheno == 0)
  my.maf <- mean(my.alleles, na.rm = TRUE)/2
  my.test2 <- glm (data = my.frame, formula = my.formula, family = binomial)
  my.p.SNP <- coef(summary(my.test2))['count.rare.variants', 4]
  
  res.data$pval.cond[ col ] <- my.p.SNP
  res.data$MAF[ col ] <- my.maf
  
  my.name <- dimnames(genotypes.inter.ichip2)[[2]][ col ]
  cat(base.code, '\t', pheno.code, '\t', sum(pheno.values == 0), '\t', sum(pheno.values == 1), '\t', my.count, '\t', covar.string, '\t', covar.p.string, '\t', my.maf, '\t', my.name, '\t', signif(my.p.SNP, 3), '\n', file = output.file.single, append = TRUE)
  cat(base.code, '\t', pheno.code, '\t', sum(pheno.values == 0), '\t', sum(pheno.values == 1), '\t', my.count, '\t', covar.string, '\t', covar.p.string, '\t', my.maf, '\t', my.name, '\t', signif(my.p.SNP, 3), '\n')
}

message("Output in ", output.file.single)

