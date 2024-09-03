dataset <- 'TCGA'
infile <- 'tcga_beta_mixture_model.rds'
# dataset <- 'Lund'
# infile <- 'lund_beta_mixture_model.rds'

indir <- '/Users/danielmonyak/Library/CloudStorage/Box-Box/PROJECT 06023: MolClocks/MolClock_Paper_1/1. Analytic Datasets/Beta Peak Decomposition'

BetaMixture <- readRDS(file.path(indir, dataset, infile))
outdir <- file.path(indir, dataset, 'readable_tables')

write.table(BetaMixture$alpha, file.path(outdir, 'alpha.txt'), quote=FALSE)
write.table(BetaMixture$delta, file.path(outdir, 'delta.txt'), quote=FALSE)
write.table(BetaMixture$phi, file.path(outdir, 'phi.txt'), quote=FALSE)

tumors <- names(BetaMixture$Z)

# maxProbGroup <- list()
# tau <- list()
# dummy.counts <- as.data.frame(matrix(rep(0, 3), nrow=1))
# names(dummy.counts) <- c('1', '2', '3')

for (tum in tumors) {
  #tum <- tumors[1]
  maxProbGroup[[tum]] <- as.data.frame(apply(BetaMixture$Z[[tum]], 1, which.max))
  names(maxProbGroup[[tum]]) <- tum

  # group.counts <- as.data.frame(t(as.matrix(table(maxProbGroup[[tum]]))))
  # group.counts.fixed <- colSums(dplyr::bind_rows(group.counts, dummy.counts), na.rm=TRUE)
  # tau[[tum]] <- group.counts.fixed / nrow(maxProbGroup[[tum]])
}
# tau.df <- t(data.frame(tau))
# write.table(tau.df, file.path(outdir, 'tau.txt'), quote=FALSE)

modes <- (BetaMixture$alpha - 1)/(BetaMixture$alpha + BetaMixture$delta - 2)
write.table(modes, file.path(outdir, 'modes.txt'), quote = F)
