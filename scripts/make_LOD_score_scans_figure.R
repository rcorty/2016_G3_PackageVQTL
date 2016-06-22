library(vqtl)

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


a1 <- scanonevar(cross = my.cross,
								 mean.formula = 'phenotype1 ~ sex + mean.QTL.add + mean.QTL.dom',
								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
								 return.models = TRUE)$varscan
a2 <- scanone(cross = my.cross,
							pheno.col = 'phenotype1')

b1 <- scanonevar(cross = my.cross,
								 mean.formula = 'phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom',
								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
								 return.models = TRUE)$varscan
b2 <- scanone(cross = my.cross,
							pheno.col = 'phenotype2')

c1 <- scanonevar(cross = my.cross,
								 mean.formula = 'phenotype3 ~ sex + mean.QTL.add + mean.QTL.dom',
								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
								 return.models = TRUE)$varscan
c2 <- scanone(cross = my.cross,
							pheno.col = 'phenotype3')

d1 <- scanonevar(cross = my.cross,
								 mean.formula = 'phenotype4 ~ sex + mean.QTL.add + mean.QTL.dom',
								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
								 return.models = TRUE)$varscan
d2 <- scanone(cross = my.cross,
							pheno.col = 'phenotype4')


pdf(file = '../2016_G3_PackageVQTL_Corty/images/LOD_scans.pdf', width = 6, height = 8)
par(mfrow = c(4, 1), mar = c(3.1, 3.1, 3.1, 2.1))
plot(x = a1, y = a2, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = c(0, 17), title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
plot(x = b1, y = b2, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = c(0, 17), title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
plot(x = c1, y = c2, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = c(0, 17), title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
plot(x = d1, y = d2, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = c(0, 17), title.cex = 1.2, line.width = 1.5)
dev.off()
