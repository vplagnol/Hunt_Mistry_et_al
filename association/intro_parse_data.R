library(snpStats)

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}



## read arguments
myArgs <- getArgs()

if ('code' %in% names(myArgs)) {code <- myArgs[[ 'code' ]];}
if ('pheno' %in% names(myArgs)) {
  pheno <- myArgs[['pheno']]
  choice.phenos <- strsplit(pheno, split = '_')[[1]]
}

pheno.code <- paste(choice.phenos, collapse = '_')
nphenos <- length(choice.phenos)
base.code <- basename(code)

### load the data
input.file <- paste(code, '.RData', sep = '')
load(input.file)

if (pheno.code == 'control_psoriasis') ichip <- pso.gwas  ##replace ichip with gwas if we deal with psoriasis samples

pheno.not.control <- subset(choice.phenos, choice.phenos != 'control')
pheno.files <- paste('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/phenotypes/', pheno.not.control, '.phe', sep = '')

if (sum(!file.exists(pheno.files)) > 0) {
  print(pheno.files)
  stop('Some files do not exist')
}

if (length(pheno.files) > 1) stop("Only one non control phenotype for now")
