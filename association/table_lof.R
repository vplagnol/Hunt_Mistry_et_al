data <- read.table('support/target_regions.tab', col.names = c('chromosome', 'start', 'end', 'code', 'gene'), stringsAsFactors = FALSE)
data <- subset(data, grepl(pattern = 'full', data$code))
data <- subset(data, ! code %in% 'YDJC_full')

#data <- subset(data,  code %in% 'UBASH3A_full')

lof.code <- c('stopgain SNV', 'frameshift deletion', 'frameshift insertion')
my.phenos <- c('t1d', 'coeliac', 'ms', 'graves', 'crohns', 'psoriasis')

my.tab <- data.frame (code = data$code, n.controls = NA,controls = NA)
for (pheno in my.phenos) {
  my.tab[, pheno] <- NA
  my.tab[, paste('n', pheno, sep = '.')] <- NA
}

row.names(my.tab) <- my.tab$code


my.categories <- c("frameshift deletion","frameshift insertion","nonframeshift deletion","nonframeshift insertion","nonsynonymous SNV","splicing SNV","stopgain SNV","synonymous SNV")
my.tab.variants <- data.frame(categories = my.categories,
                              overall = 0,
                              rare = 0,
                              novel = 0)

my.freq.variants <- data.frame(categories = c('singleton', 'doubleton'),
                              overall = 0,
                              rare = 0,
                              novel = 0)


my.type.variants <- data.frame(categories = c('transition', 'transversion'),
                              overall = 0,
                              rare = 0,
                              novel = 0)


all.coding <- data.frame()
all.candidates <- data.frame()
all.rare <- data.frame()
all.novel <- data.frame()



multi.alleles <- rep(0, 4)


for (code in data$code) {
  message(code)

  input.file <- paste('data/genotypes/', code, '.RData', sep = '')
  load(input.file)

  multi.alleles <- multi.alleles + table(factor(summ$nalleles, levels = 1:4))
  
  my.tab.variants$overall <- my.tab.variants$overall +  as.numeric(table(factor(summ$ExonicFunc, levels = my.categories)))
  my.tab.variants$rare <- my.tab.variants$rare +  as.numeric(table(factor(subset(summ$ExonicFunc, summ$rare), levels = my.categories)))
  my.tab.variants$novel <- my.tab.variants$novel +  as.numeric(table(factor(subset(summ$ExonicFunc, summ$novel), levels = my.categories)))
  
  my.freq.variants[1, 'overall'] <- my.freq.variants[1, 'overall'] + sum(summ$nhets == 1 & !is.na(summ$ExonicFunc), na.rm = TRUE)
  my.freq.variants[2, 'overall'] <-  my.freq.variants[2, 'overall'] + sum(summ$nhets == 2 & !is.na(summ$ExonicFunc),  na.rm = TRUE)
  my.freq.variants[1, 'rare'] <- my.freq.variants[1, 'rare'] + sum(summ$nhets == 1 & summ$rare & !is.na(summ$ExonicFunc), na.rm = TRUE)
  my.freq.variants[2, 'rare'] <-  my.freq.variants[2, 'rare'] + sum(summ$nhets == 2 & summ$rare  & !is.na(summ$ExonicFunc), na.rm = TRUE)
  my.freq.variants[1, 'novel'] <- my.freq.variants[1, 'novel'] + sum(summ$nhets == 1 & summ$novel  & !is.na(summ$ExonicFunc), na.rm = TRUE)
  my.freq.variants[2, 'novel'] <-  my.freq.variants[2, 'novel'] + sum(summ$nhets == 2 & summ$novel  & !is.na(summ$ExonicFunc), na.rm = TRUE)


  my.mut <- paste(summ$Ref, summ$Alt, sep = '')
  ti.tv <- rep(NA, length(my.mut))
  my.mut <- paste(summ$REF, summ$ALT, sep = '')
  ti.tv <- ifelse (my.mut %in% c('AG', 'GA', 'CT', 'TC'), 1, ti.tv)
  ti.tv <- ifelse (my.mut %in% c('AC', 'CA', 'AT', 'TA', 'GT', 'TG', 'CG', 'GC'), 2, ti.tv)

  all.coding <- rbind.data.frame(all.coding, summ[, c('chromosome', 'start', 'REF', 'ALT', 'Gene')])
  all.candidates <- rbind.data.frame(all.candidates, subset(summ[, c('chromosome', 'start', 'REF', 'ALT', 'Gene')], summ$candidate))
  all.rare <- rbind.data.frame(all.rare, subset(summ[, c('chromosome', 'start', 'REF', 'ALT', 'Gene')], summ$rare))
  all.novel <- rbind.data.frame(all.novel, subset(summ[, c('chromosome', 'start', 'REF', 'ALT', 'Gene')], summ$novel))

  summ$weights <- 1
  #summ$weights <- ifelse (summ$nalleles == 1, 1, 0)
  #summ$weights <- ifelse (summ$nhets > 1 , 1, 0)

  my.type.variants[1, 'overall'] <- my.type.variants[1, 'overall'] + sum(summ$weights*(ti.tv == 1 & !is.na(summ$ExonicFunc)), na.rm = TRUE)
  my.type.variants[2, 'overall'] <-  my.type.variants[2, 'overall'] + sum(summ$weights*(ti.tv == 2 & !is.na(summ$ExonicFunc)), na.rm = TRUE)
  my.type.variants[1, 'rare'] <- my.type.variants[1, 'rare'] + sum(summ$weights*(ti.tv == 1 & summ$rare & !is.na(summ$ExonicFunc)), na.rm = TRUE)
  my.type.variants[2, 'rare'] <-  my.type.variants[2, 'rare'] + sum(summ$weights*(ti.tv == 2 & summ$rare  & !is.na(summ$ExonicFunc)), na.rm = TRUE)
  my.type.variants[1, 'novel'] <- my.type.variants[1, 'novel'] + sum(summ$weights*(ti.tv == 1 & summ$novel  & !is.na(summ$ExonicFunc)), na.rm = TRUE)
  my.type.variants[2, 'novel'] <-  my.type.variants[2, 'novel'] + sum(summ$weights*(ti.tv == 2 & summ$novel  & !is.na(summ$ExonicFunc)), na.rm = TRUE)

  message(sum(ti.tv == 1 & summ$novel  & !is.na(summ$ExonicFunc), na.rm = TRUE) /  sum(ti.tv == 2 & summ$novel  & !is.na(summ$ExonicFunc), na.rm = TRUE), '   ', sum(ti.tv == 1 & summ$novel  & !is.na(summ$ExonicFunc), na.rm = TRUE))
  
  n.lof <- sum(summ$ExonicFunc %in% lof.code > 0, na.rm = TRUE)

  if (n.lof > 0) {
    for (pheno in my.phenos) {
      pheno.files <- paste('/data_n2/hmw208/Fluidigm_resequencing/1014AAP11O1/libs_01to34_pseq/phenotypes/', pheno, '.phe', sep = '')
      pheno.table <- read.table(pheno.files, stringsAsFactors = FALSE)
      
      if (pheno == my.phenos[1]) { ##the first time I use controls
        controls <- subset(pheno.table, V2 == 1, 'V1', drop = TRUE)
        geno.loc.ctl <- genotypes[ summ$ExonicFunc %in% lof.code, dimnames(genotypes)[[2]] %in% controls, drop = FALSE ]
        my.tab[code, 'controls' ] <- sum(geno.loc.ctl, na.rm = TRUE)
        my.tab[code, 'n.controls' ] <- ncol(geno.loc.ctl)
      }
      
      cases <- subset(pheno.table, V2 == 2, 'V1', drop = TRUE)      
      geno.loc <- genotypes[ summ$ExonicFunc %in% lof.code, dimnames(genotypes)[[2]] %in% cases, drop = FALSE ]
      my.tab[code, pheno ] <- sum(geno.loc, na.rm = TRUE)
      my.tab[code, paste('n', pheno, sep = '.') ] <- ncol(geno.loc)


      my.n <- c(ncol(geno.loc.ctl), ncol(geno.loc))
      my.counts <- c(  my.tab[code, 'controls' ] ,   my.tab[code, pheno ] )
        
      if (sum(my.counts) > 0) {
        my.p <- fisher.test( matrix (data = c(my.n, my.counts), ncol = 2, nrow = 2))$p.value
        if (my.p < 0.05) {
          message(code, ' ', pheno, ' ', my.p)
        }
      }
      
      
    }
  } else {
    for (pheno in c('controls', my.phenos)) {my.tab[ code, pheno ] <- 0}
   }
}


pretty.tab <- my.tab[,c('controls', my.phenos)]
my.ns <- paste('(n=', as(my.tab[1,paste('n', c('controls', my.phenos), sep = '.')], 'numeric'), ')', sep = '')
names(pretty.tab) <- paste(names(pretty.tab), my.ns)
row.names(pretty.tab) <- gsub(pattern = '_full', replacement = '', row.names(pretty.tab) )

write.csv(x = pretty.tab, file = 'LOF_variants.csv', row.names = TRUE)



#######
fram.indels <- subset(my.tab.variants, categories ==  'frameshift deletion') +  subset(my.tab.variants, categories ==  'frameshift insertion') 
fram.indels[[1]] <- 'frameshift indels'


nfram.indels <- subset(my.tab.variants, categories ==  'nonframeshift deletion') +  subset(my.tab.variants, categories ==  'nonframeshift insertion') 
nfram.indels[[1]] <- 'nonframeshift indels'

my.tab.variants <- rbind.data.frame(my.tab.variants, rbind.data.frame(fram.indels, nfram.indels))
my.tab.variants <- subset(my.tab.variants, categories %in% c('frameshift indels', 'nonframeshift indels', 'synonymous SNV', 'nonsynonymous SNV', 'splicing SNV', 'stopgain SNV'))

my.ti.tv <- my.type.variants[1,] /my.type.variants[2,]
my.ti.tv[[1]] <- 'transition transversion ratio'

complete.table <- rbind.data.frame(my.tab.variants, my.freq.variants)
complete.table[1,] <- as(complete.table[,1], 'character')
complete.table <- rbind.data.frame(my.tab.variants, my.ti.tv)

write.csv(x = complete.table, file = 'count_variants.csv', row.names = TRUE)


write.table(x = all.coding, row.names = FALSE, file = 'all_coding.tab', sep = '\t', col.names = FALSE, quote = FALSE)
write.table(x = all.candidates, row.names = FALSE, file = 'all_candidate.tab', sep = '\t', col.names = FALSE, quote = FALSE)
write.table(x = all.rare, row.names = FALSE, file = 'all_rare.tab', sep = '\t', col.names = FALSE, quote = FALSE)
write.table(x = all.novel, row.names = FALSE, file = 'all_novel.tab', sep = '\t', col.names = FALSE, quote = FALSE)
