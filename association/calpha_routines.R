
get.count.table <- function(x, pheno, signatures, combine.singletons = TRUE) {

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
  if (combine.singletons) table.results <- rbind(table.results, singletons)
  return(table.results)
}


calpha.test <- function(x, proportions, phenos = c('control', 'coeliac'), verbose = FALSE ){

  x <- x[, phenos]
  proportions <- proportions[ phenos ]
  total.counts <- apply(x, FUN = sum, MAR = 1)
  U.loc <- (x[,1] - total.counts*proportions[1])^2 - total.counts*proportions[1]*(1-proportions[1])
  U <- sum(U.loc) 
  return (U)
}





calpha.test.multi <- function(x, proportions, phenos = c('control', 'coeliac'), verbose = FALSE ){

  nphenos <- length(phenos)
  x <- x[, phenos]
  proportions <- proportions[ phenos ]
  U <- rep(0, times = nphenos - 1)

  for (i in 1:(nphenos - 1)) {
### first clean up the matrix of counts
    x.loc <- x
   
    total.counts <- apply(x.loc, FUN = sum, MAR = 1)
    new.single <- total.counts == 1 & dimnames(x.loc)[[1]] != 'singletons'
    if (sum(new.single) > 0) {x.loc['singletons',] <- x.loc['singletons',] + apply(FUN = sum, x.loc[ new.single, , drop = FALSE], MAR = 2)}

    x.loc <- x.loc[ total.counts >= 2,]
    total.counts <- apply(x.loc, FUN = sum,MAR = 1)

    ### adjust to make it binomial if we have more than 2 cohorts
    other.cohorts <-  apply(MAR = 1, FUN = sum, x.loc[, -1, drop = FALSE])
    x.loc <- x.loc[, 1, drop = FALSE]
    x.loc <- cbind(x.loc, other.cohorts)
    proportions.loc <- proportions[ phenos]
    proportions.loc <- c(proportions.loc[ 1 ], sum(proportions.loc[ -1 ]))
    proportions.loc <- proportions.loc/sum(proportions.loc)

    ### compute the expected binomial proportions     
    U.loc <- (x.loc[,1] - total.counts*proportions.loc[1])^2 - total.counts*proportions.loc[1]*(1-proportions.loc[1])
    U[i] <- U[i] + sum(U.loc) 
    
    if (verbose) {
      print(proportions.loc)
      print(cbind(x.loc, U.loc))
    }
    
    ### Now remove the phenotype we just tested
    phenos <- phenos[-1]
    x <- x[, -1, drop = FALSE]
    proportions <- proportions[-1]
  }

  return(U)
}



