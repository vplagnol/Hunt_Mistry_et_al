

data <- read.table('depth.tab')
mat <- data.frame()
for (d in c(20:250)) {
  mat <- rbind.data.frame(mat, c(d, sum(subset(data, V2 %in% as.character(d:(d+250)), 'V1', drop = TRUE))/sum(data$V1)) )
}

total.nb.calls <- 41911*4577

pdf('depth_vs_proportion.pdf', width = 5, height = 5.5)

plot (x = mat[,1],
      y = mat[,2],
      xlab = 'Read depth',
      ylab = 'Proportion of non-missing genotype calls with depth > x',
      type = 'l')


dev.off()

print(sum(as.numeric(as.character(data$V1)), na.rm = TRUE))
