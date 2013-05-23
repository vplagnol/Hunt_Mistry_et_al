library(snpStats)
library(VPlib)

getArgs <- function() {
  myargs.list <- strsplit(grep("=",gsub("--","",commandArgs()),value=TRUE),"=")
  myargs <- lapply(myargs.list,function(x) x[2] )
  names(myargs) <- lapply(myargs.list,function(x) x[1])
  return (myargs)
}

checks <- TRUE

my.func <- c('frameshift insertion', 'frameshift deletion', 'nonframeshift deletion',  'nonframeshift insertion', 'nonsynonymous SNV', 'stopgain SNV', 'stoploss SNV', 'splicing SNV')
my.lof.func <- c('frameshift insertion', 'frameshift deletion', 'stopgain SNV')


code <- 'data/genotypes/ZFP36L1_full'
myArgs <- getArgs()
if ('code' %in% names(myArgs)) {code <- myArgs[[ 'code' ]]}

loc.gene <- gsub(basename(code), pattern = '_full', replacement = '')

sampleIDs.file <- paste(code, '_IDs.tab', sep = '')

IDs <- scan(sampleIDs.file, what = 'character')



################################## parse annovar file

my.names <- c('Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','ESP6500_ALL','X1000g2012apr_ALL','dbSNP137','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LJB_MutationTaster','LJB_MutationTaster_Pred','LJB_GERP++','Chr','Start','End','Ref','Obs','Call', 'junk')

#my.names <- c('Func','Gene','ExonicFunc','AAChange','Conserved','SegDup','1000g2010nov_ALL','1000g2012feb_ALL','1000g2012apr_ALL','dbSNP130','AVSIFT','LJB_PhyloP','LJB_PhyloP_Pred','LJB_SIFT','LJB_SIFT_Pred','LJB_PolyPhen2','LJB_PolyPhen2_Pred','LJB_LRT','LJB_LRT_Pred','LRT_MutationTaster','LRT_MutationTaster_Pred','LJB_GERP++','ESP6500_ALL','dbSNP135','dbSNP135Flagged','dbSNP135Common','dbSNP135Mult','cg46','cg69','ESP6500_EA','ESP6500_AA','Omim','Chr','Start','End','Ref','Obs','Call', 'QUAL', 'Depth')

annovar <- read.csv('support/fluidigm50k_final15jan2013_noIFCbarcode__has_phenotype.indiv_4377sites_41911indiv.genome_summary.csv', col.names = my.names, stringsAsFactor = FALSE, na.string = c('NA', ''))
annovar$Chr <- gsub(annovar$Chr, pattern = 'chr', replacement = '') 
annovar$signature <- paste(annovar$Chr, annovar$Start, annovar$Ref, annovar$Obs, sep = '_')

##better handle of splice variants
annovar$ExonicFunc <-  ifelse ( grepl(annovar$Func, pattern = 'splicing') & annovar$ExonicFunc %in% c(NA, 'nonsynonymous SNV', 'synonymous SNV'), 'splicing SNV', annovar$ExonicFunc)


annovar$rare <- (is.na(annovar$X1000g2012apr_ALL) | annovar$X1000g2012apr_ALL < 0.005) &  (is.na(annovar$ESP6500_ALL) | annovar$ESP6500_ALL < 0.005)
annovar$novel <- is.na(annovar$X1000g2012apr_ALL) & is.na(annovar$ESP6500_ALL) & is.na(annovar$dbSNP137)
annovar$candidate <- annovar$rare & (  (annovar$ExonicFunc %in% my.func) | (annovar$Func %in% c('splicing', 'exonic;splicing')) )


just.once <- FALSE
if (just.once) {
  coding <- subset(annovar, ! is.na(annovar$ExonicFunc))
  write.table(x = coding[, c('Chr', 'Start', 'Ref', 'Obs')], row.names = FALSE, quote = FALSE, sep = ' ', file = 'all_coding_variants.tab')

  rare <-  subset(coding, rare)
  write.table(x = coding[, c('Chr', 'Start', 'Ref', 'Obs')], row.names = FALSE, quote = FALSE, sep = ' ', file = 'all_rare_variants.tab')
  
  candidate <- subset(annovar, candidate)
  write.table(x = candidate[, c('Chr', 'Start', 'Ref', 'Obs')], row.names = FALSE, quote = FALSE, sep = ' ', file = 'all_candidate_variants.tab')
  stop()
}



#write.csv(x = subset(annovar, candidate), file = 'annovar_candidates.csv', row.names = FALSE); stop()



summ <- read.table(paste(code, '_summary.tab', sep = ''), header = TRUE, stringsAsFactors = FALSE)

nSNPs <- nrow(summ)
nsamples <- length(IDs)

summ$signature <- paste(summ$chromosome, summ$start, summ$REF, summ$ALT, sep = '_')

###################################### fix the signature for the variants
for (spos in c('A', 'C', 'G', 'T')) {
  to.fix <- grepl( pattern = paste('^', spos, sep = ''), x = summ$REF) &  grepl( pattern = paste('^', spos, sep = ''), x = summ$ALT)
  summ$REF <- ifelse (to.fix, gsub(pattern = paste('^', spos, sep = ''), replacement = '', x = summ$REF), summ$REF)
  summ$ALT <- ifelse (to.fix, gsub(pattern = paste('^', spos, sep = ''), replacement = '', x = summ$ALT), summ$ALT)
}


summ$REF <- ifelse (summ$REF == '', '-', summ$REF)
summ$ALT <- ifelse (summ$ALT == '', '-', summ$ALT)


summ$signature <- paste(summ$chromosome,
                                      ifelse(summ$ALT == '-', summ$start + 1, summ$start),
                                      summ$REF,
                                      summ$ALT,
                                      sep = '_')


########### extract the region from the iChip files
chromosome <- summ$chromosome[1]
my.range.kb <- round(   (range(summ$start) + c(-200000, + 200000))/1000 )

ichip.files <- '/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/fluidigm50k_final15jan2013_noIFCbarcode.i-view_indiv_withphenotype_plink__toload'

my.cmd<- paste('/data_n2/vplagnol/Software/plink-1.07-x86_64/plink --noweb --chr ', chromosome, '--from-kb', my.range.kb[1], '--to-kb',  my.range.kb[2], '--bfile', ichip.files, '--make-bed --out ', code)
system(my.cmd)
ichip <- read.plink(code)

###########
pso.gwas <- read.plink('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/psoriasis_CCC2_GWAS/CCC2_GWAS_psoriasis_3')
unmapped.pos <- read.table('data/pso_gwas_remove.tab')$V2
good <- ! pso.gwas$map$position %in% unmapped.pos
pso.gwas$map <- pso.gwas$map[good,]
pso.gwas$genotypes <- pso.gwas$genotypes[, good ]

pso.gwas$map$position <- hg18.to.hg19 (chrom = pso.gwas$map$chromosome,
                               position = pso.gwas$map$position,
                               liftOver = "/data_n2/vplagnol/Data/UCSC/liftOver",
                               chain = "/data_n2/vplagnol/Data/UCSC/hg18ToHg19.over.chain.gz")

n.good <- pso.gwas$map$chromosome == chromosome & pso.gwas$map$position > my.range.kb[1]*1000 &  pso.gwas$map$position < my.range.kb[2]*1000
pso.gwas$map <- pso.gwas$map[n.good,]
pso.gwas$genotypes <- pso.gwas$genotypes[, n.good ]


####################### Now read the genotype data
my.data <- scan( paste(code, '_geno.tab', sep = ''), sep = '\t')
genotypes <- matrix(data = my.data, ncol = nsamples, nrow = nSNPs, byrow = TRUE)
dimnames(genotypes) <- list(summ$signature, IDs)

summ$prop.missing <- apply(MAR = 1, FUN = function(x) {sum(is.na(x))}, genotypes) / ncol(genotypes)



##################################################################################### merge with annovar and define candidates
print(table(summ$signature  %in% annovar$signature))
print(subset(summ$signature, ! summ$signature  %in% annovar$signature))

summ <- merge(summ, annovar[, c('Func','Gene','ExonicFunc','AAChange', 'X1000g2012apr_ALL', 'dbSNP137', 'ESP6500_ALL', 'signature', 'candidate', 'rare', 'novel')], by = 'signature', all.x = TRUE, sort = FALSE)
summ <- summ[ match(table = summ$signature, dimnames(genotypes)[[1]]), ]  ##reorder properly
if (sum(dimnames(genotypes)[[1]] != summ$signature) > 0) stop()



###### remove non coding variants
keep.variants <- !is.na(summ$ExonicFunc) & grepl(pattern = loc.gene, summ$Gene)
summ <- summ[ keep.variants, ]
genotypes <- genotypes[ keep.variants,]

output.file.summary <- paste('data/gene_summary/', loc.gene, '.tab', sep  = '')
cat(code, '\t', loc.gene, '\t', nrow(summ), '\t', sum(summ$candidate), '\n', file = output.file.summary)

### Now look for problematic SNPs
problem <- subset(summ, candidate & MAF.nic > 0.02)
if (nrow(problem) > 0) {
  output.problem <- paste(code, '.problem', sep = '')
  write.table(x = problem, sep = '\t', row.names = FALSE, quote = FALSE, file = output.problem)
}




######## Now run some checks
if (checks) {
  print(table(dimnames(ichip$genotypes)[[1]] %in% IDs)) #should be all TRUE
  
  high.freq <- subset(summ, MAF > 0.3) ## pick the common SNPs
  if (nrow(high.freq) > 0) {
    
    for (j in 1:nrow(high.freq)) {
      my.pos <- high.freq$start[ j ]
      pos.ichip <- which (ichip$map$position == my.pos)
      ichip.geno <- as(ichip$genotypes[, pos.ichip], 'numeric')[ match(IDs, table = dimnames(ichip$genotypes)[[1]]) ]	
      print(table(ichip.geno, genotypes[ high.freq$signature[j],]))
      print(cor(ichip.geno, genotypes[ high.freq$signature[j],], use = 'pairwise.complete.obs'))
    }
  }
}


####### save the data
output.file <- paste(code, '.RData', sep = '')
save(list = c('pso.gwas', 'ichip', 'genotypes', 'nSNPs', 'nsamples', 'summ'), file = output.file)


