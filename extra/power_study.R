
n.controls <- 18000
n.cases <- 2100
odds.ratio <- 4

power.rare <- function(MAF.controls, odds.ratio, n.controls, n.cases, pvalue = 10^(-4)) {

  n.rare.controls <- 2*n.controls*MAF.controls
  n.rare.cases <- 2*n.cases*odds.ratio*MAF.controls

  my.sd<- sqrt(sum(c(1/n.rare.controls, 1/n.rare.cases, 1/(2*n.controls - n.rare.controls), 1/(2*n.cases - n.rare.cases))))


  ncpv <- (log(odds.ratio)/my.sd)^2
  
  threshold <- qchisq(p = 1 - pvalue, df = 1)
  return ( pchisq(q = threshold, df = 1, ncp = ncpv, lower.tail = FALSE))
}







output.pdf <- 'power_study.pdf'
pdf(output.pdf, width = 8, height = 7)
par(mfrow = c(2, 2))

for (f.controls in c(0.001, 0.005, 0.01, 0.05)) {
  
  if (f.controls %in% c(0.001)) my.ORs <- seq(1., 5, by = 0.01)
  if (f.controls %in% c(0.005)) my.ORs <- seq(1., 3, by = 0.01)
  if (f.controls %in% c(0.01)) my.ORs <- seq(1., 2, by = 0.01)
  if (f.controls %in% c(0.05)) my.ORs <- seq(1., 1.6, by = 0.01)
  
  my.power.2100 <- rep(NA, length(my.ORs))
  my.power.6500 <- rep(NA, length(my.ORs))
  
  for (i in 1:length(my.ORs)) {
    my.power.2100[i] <- power.rare (MAF.controls = f.controls, odds.ratio = my.ORs[i], n.cases = 2100, n.controls = 17000, pvalue = 10^-4)
    my.power.6500[i] <- power.rare (MAF.controls = f.controls, odds.ratio = my.ORs[i], n.cases = 6500, n.controls = 17000, pvalue = 10^-4)
  }
  
  plot (x = my.ORs,
        y = my.power.6500,
        col = 'black',
        ylim = c(0,1),
        type = 'l',
        xlab = 'Odds ratio',
        ylab = expression(paste('Probability p <', 10^-4)),
        main = paste('Allele frequency in controls ', f.controls*100, '%', sep = ''),
        sub = '17K controls and 6.5K/2.1K cases (black/red)')
  
  lines(x = my.ORs,
        y = my.power.2100,
        col = 'red',
        type = 'l')
  
  abline(h = 0.9)
  
  abline(v = my.ORs[ which.min(abs(my.power.2100 - 0.9)) ], col = 'red')
  abline(v = my.ORs[ which.min(abs(my.power.6500 - 0.9)) ], col = 'black')
}
dev.off()


print(output.pdf)
