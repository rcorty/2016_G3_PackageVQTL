library(tidyr); library(dplyr); library(vqtl)

set.seed(27599)

my.cross <- sim.cross(map = sim.map(len = rep(100, 4), n.mar = 30, eq.spacing = TRUE, include.x = FALSE),
											n.ind = 200,
											type = 'f2')
my.cross$pheno$sex <- rep(x = c(0, 1), each = 100)
my.cross <- calc.genoprob(my.cross)

my.cross$pheno$phenotype1 <- rnorm(n = nind(my.cross))
my.cross$pheno$phenotype2 <- rnorm(n = nind(my.cross), mean = 0.8*my.cross$geno$`1`$data[,15])
my.cross$pheno$phenotype3 <- rnorm(n = nind(my.cross), sd = my.cross$geno$`2`$data[,15])
my.cross$pheno$phenotype4 <- rnorm(n = nind(my.cross), mean = my.cross$geno$`3`$data[,15], sd = my.cross$geno$`3`$data[,15])

