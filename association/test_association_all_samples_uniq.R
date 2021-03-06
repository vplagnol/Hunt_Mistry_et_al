

##default parameters
choice.phenos <- c('control', 'coeliac')
code <- 'data/genotypes/NCF2_full'
#code <- 'data/genotypes/IL18RAP_full'

source('scripts/association/intro_parse_data.R')

get.count.table.uniq <- function(x, pheno, signatures, combine.singletons = TRUE) {

  n.phenos <- length(levels(pheno))
  table.results <- matrix(nrow = 0, ncol = n.phenos)
  singletons <- rep(0, n.phenos)

  for (i in 1:nrow(x)) {  #for each variant in that gene
    counts <- tapply(x[i,], IND = pheno, FUN = sum, na.rm = TRUE)
    total.counts <- sum(counts)
    if (combine.singletons && (total.counts == 1)) {singletons <- singletons + counts} else {
      table.results <- rbind(table.results, counts)
      dimnames(table.results)[[1]][ nrow(table.results) ] <- signatures[ i ]
    }
  }
  table.results <- rbind(table.results, singletons)
  return(table.results)
}


#### define output files
output.folder.UNIQ <- 'data/pvalues_all_samples_UNIQ'
output.file.UNIQ <- paste(output.folder.UNIQ, '/', base.code, '_', pheno.code, '.out', sep = '')
output.table.UNIQ <- paste(output.folder.UNIQ, '/', base.code, '_', pheno.code, '.tab', sep = '')
output.sim.UNIQ <- paste(output.folder.UNIQ, '/', base.code, '_', pheno.code, '.RData', sep = '')

if (!file.exists(output.folder.UNIQ )) dir.create(output.folder.UNIQ)


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

#diallelic <- summ$nalleles == 1 & summ$candidate
diallelic <- summ$candidate
geno.di <- t(genotypes.full[ , diallelic ])
summ.di <- summ[ diallelic, ]

pheno.full.loc <- factor( choice.phenos[ my.frame$pheno + 1 ], levels = choice.phenos)
proportions <-  table(pheno.full.loc)/length(pheno.full.loc)

source('scripts/association/calpha_routines.R')
my.counts <- get.count.table (geno.di, pheno = pheno.full.loc, summ.di$signature, combine.singletons = FALSE)


my.T.cases <- sum(geno.di[apply(geno.di[, pheno.full.loc == 'control'], MAR = 1, FUN = sum, na.rm = TRUE) == 0,], na.rm = TRUE)
my.T.controls <- sum(geno.di[apply(geno.di[, pheno.full.loc == choice.phenos[2] ], MAR = 1, FUN = sum, na.rm = TRUE) == 0,], na.rm = TRUE)
message('Uniq count: ', my.T.cases, '\t', my.T.controls)

##### Now simulate extensively
nsim <- 10000
my.T.sim.cases <- c()
my.T.sim.controls <- c()

for (j in 1:nsim) {
  set.seed(j)
  new.pheno <-  sample(pheno.full.loc, replace = FALSE)
  not.in.controls <- sum(geno.di[apply(geno.di[, new.pheno == 'control'], MAR = 1, FUN = sum, na.rm = TRUE) == 0,], na.rm = TRUE)
  not.in.cases <- sum(geno.di[apply(geno.di[, new.pheno == choice.phenos[2]], MAR = 1, FUN = sum, na.rm = TRUE) == 0,], na.rm = TRUE)
  my.T.sim.cases <- c(my.T.sim.cases, not.in.controls)
  my.T.sim.controls <- c(my.T.sim.controls, not.in.cases)
  message(j, ' ', my.T.sim.cases[j], '\t', my.T.sim.controls[j])
}

my.p.safe.cases <- (0.5*sum(my.T.sim.cases == my.T.cases) + sum( my.T.sim.cases > my.T.cases))/nsim
my.p.safe.controls <- (0.5*sum(my.T.sim.controls == my.T.controls) + sum( my.T.sim.controls > my.T.controls))/nsim


if (my.p.safe.cases == 0) my.p.safe.cases <- 0.5/nsim
if (my.p.safe.controls == 0) my.p.safe.controls <- 0.5/nsim

save(list = c('my.T.sim.cases', 'my.T.sim.controls', 'my.T.cases', 'my.T.controls'), file = output.sim.UNIQ)

       
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', nsim, '\t', signif(my.p.safe.cases, 3), '\t', signif(my.p.safe.controls, 3), '\n', file = output.file.UNIQ)
cat(base.code, '\t', pheno.code, '\t', sum(my.frame$pheno == 0), '\t', sum(my.frame$pheno == 1), '\t', nsim, '\t', signif(my.p.safe.cases, 3), '\t', signif(my.p.safe.controls, 3), '\n')

