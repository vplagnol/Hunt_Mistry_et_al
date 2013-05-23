##default parameters
choice.phenos <- c('control', 't1d')
code <- 'data/genotypes/IL2_full'

source('scripts/association/intro_parse_data.R')

output.file <- paste('data/combined_stepwise/', base.code, '_', pheno.code, '.out', sep = '')
output.details <- paste('data/combined_stepwise/', base.code, '_', pheno.code, '.details', sep = '')


###### merge with ichip
good.snps <- (summ$nalleles == 1) &  ! (summ$start %in% ichip$map$position)
summ.restricted <- summ[ good.snps, ] ## this is the fluidigm support
genotypes <- t(genotypes[good.snps, dimnames(ichip$genotypes)[[1]] ])
dimnames(genotypes)[[2]] <- paste('chr', dimnames(genotypes)[[2]], sep = '')
fluidigm <- new('SnpMatrix', .Data = genotypes + 1)


combined.snpStats <- snp.cbind(x = ichip$genotypes, y = fluidigm)
combined.support <- rbind.data.frame (data.frame(SNP = ichip$map$snp.name, chromosome = ichip$map$chromosome, position = ichip$map$position),
                                      data.frame(SNP = paste('chr', summ.restricted$signature, sep = ''), chromosome = summ.restricted$chromosome, position = summ.restricted$start))
combined.support$coding <- combined.support$position %in% summ$start
row.names(combined.support) <- combined.support$SNP

####### now the phenotype data
pheno.table <- read.table(pheno.files[1])
pheno.values <- pheno.table$V2[ match(dimnames(combined.snpStats)[[1]], table = pheno.table$V1) ] -1

my.frame <- data.frame(pheno = pheno.values)
row.names(my.frame) <- dimnames(combined.snpStats)[[1]]
good.samples <- !is.na(my.frame$pheno)

combined.snpStats <- combined.snpStats[good.samples, ]
combined.num <- as(combined.snpStats, 'numeric')
#print(table(dimnames(combined.num)[[1]] == row.names(my.frame)))
my.frame <- subset(my.frame, good.samples)


### another SNP QC, to remove the super rare ones
my.counts <- apply(combined.num, FUN = sum, MAR = 2, na.rm = TRUE)
common.SNPs <- my.counts > 3

combined.snpStats <- combined.snpStats[, common.SNPs]
combined.support <- combined.support[ common.SNPs,]
my.sum <- col.summary(combined.snpStats)
combined.support$MAF <- my.sum$MAF





###### Now the stepwise regression
nSNPs <- ncol(combined.num)
min.P <- 0
nrounds <- 0

my.formula <- 'pheno ~ 1'
my.snps <- c()
my.pvals <- c()



while (min.P < 10^(-4)) {

  my.tests <- snp.rhs.tests (snp.data = combined.snpStats,
                             formula = as.formula(my.formula),
                             family = 'binomial',
                             data = my.frame)

  my.p <- p.value(my.tests)

  best.SNP <- which.min(my.p)
  best.SNP.name <- dimnames(combined.snpStats)[[2]][ best.SNP ]
  min.P <- my.p[ best.SNP ]

  
  if (min.P < 10^-4) {
    my.snps <- c(my.snps, best.SNP.name)
    my.pvals <- c(my.pvals, min.P)
    my.formula <- paste(my.formula, '+', best.SNP.name)
    my.frame[, best.SNP.name] <- as(as(combined.snpStats[, best.SNP], 'numeric'), 'numeric')
    
    message(my.formula, ' ', min.P)
  }
}

### create the file with the details
details <- data.frame(SNP = my.snps,
                      stepwise.pval = my.pvals,
                      MAF = my.sum[ my.snps, 'MAF'])

details$coding <- combined.support[ my.snps, 'coding']


is.fluidigm <- paste('fluidigm', sum(grepl(pattern = '^chr', my.snps)) > 0, sep = '')
is.coding <- paste('coding', sum(details$coding) > 0, sep = '')


##### now output the whole thing
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', is.coding, '\t', is.fluidigm, '\t', my.formula, '\n')
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', is.coding, '\t', is.fluidigm, '\t', my.formula, '\n', file = output.file)

write.table(x = details, file = output.details, row.names = FALSE, quote = FALSE, sep = '\t')
print(output.details)
