library(tidyr); library(dplyr); library(vqtl)

set.seed(27599)

my.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 20, eq.spacing = TRUE, include.x = FALSE),
													 n.ind = 400,
													 type = 'f2')
my.cross$pheno$sex <- rep(x = c(0, 1), each = 200)
my.cross <- qtl::calc.genoprob(my.cross)

my.cross$pheno$phenotype1 <- rnorm(n = qtl::nind(my.cross))
my.cross$pheno$phenotype2 <- rnorm(n = qtl::nind(my.cross), mean = 0.6*my.cross$geno$`1`$data[,10])
my.cross$pheno$phenotype3 <- rnorm(n = qtl::nind(my.cross), sd = exp(0.35*(my.cross$geno$`2`$data[,10] - 2)))
my.cross$pheno$phenotype4 <- rnorm(n = qtl::nind(my.cross), mean = 0.3*my.cross$geno$`3`$data[,10], sd = exp(0.25*(my.cross$geno$`3`$data[,10] - 2)))


# make fig 1 -- LOD score scans
ymax <- 30

a0 <- qtl::scanone(cross = my.cross,
									 pheno.col = 'phenotype1')
a1 <- scanonevar(cross = my.cross,
								 mean.formula = phenotype1 ~ sex + mean.QTL.add + mean.QTL.dom,
								 var.formula = ~sex + var.QTL.add + var.QTL.dom)
p <- plot(x = a1, y = a0, ylim = c(0, ymax))
p
ggplot2::ggsave(filename = 'images/LOD_scan_phen1.pdf', height = 3, width = 9)


b0 <- qtl::scanone(cross = my.cross,
									 pheno.col = 'phenotype2')
b1 <- scanonevar(cross = my.cross,
								 mean.formula = phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom,
								 var.formula = ~sex + var.QTL.add + var.QTL.dom)
p <- plot(x = b1, y = b0, ylim = c(0, ymax))
p
ggplot2::ggsave(filename = 'images/LOD_scan_phen2.pdf', height = 3, width = 9)


c0 <- qtl::scanone(cross = my.cross,
									 pheno.col = 'phenotype3')
c1 <- scanonevar(cross = my.cross,
								 mean.formula = phenotype3 ~ sex + mean.QTL.add + mean.QTL.dom,
								 var.formula = ~sex + var.QTL.add + var.QTL.dom)
p <- plot(x = c1, y = c0, ylim = c(0, ymax))
p
ggplot2::ggsave(filename = 'images/LOD_scan_phen3.pdf', height = 3, width = 9)


d0 <- qtl::scanone(cross = my.cross,
									 pheno.col = 'phenotype4')
d1 <- scanonevar(cross = my.cross,
								 mean.formula = phenotype4 ~ sex + mean.QTL.add + mean.QTL.dom,
								 var.formula = ~sex + var.QTL.add + var.QTL.dom)
p <- plot(x = d1, y = d0, ylim = c(0, ymax))
p
ggplot2::ggsave(filename = 'images/LOD_scan_phen4.pdf', height = 3, width = 9)




# make fig 2 -- empirircal p-value scans
n.perms <- 50
a2 <- scanonevar.perm(sov = a1, n.perms = n.perms)
b2 <- scanonevar.perm(sov = b1, n.perms = 2*n.perms)
c2 <- scanonevar.perm(sov = c1, n.perms = 2*n.perms)
d2 <- scanonevar.perm(sov = d1, n.perms = 2*n.perms)

ymax <- 10
plot(x = a2, ylim = c(0, ymax))
ggplot2::ggsave(filename = 'images/empir_p_scan_phen1.pdf', height = 3, width = 9)

plot(x = b2, ylim = c(0, ymax))
ggplot2::ggsave(filename = 'images/empir_p_scan_phen2.pdf', height = 3, width = 9)

plot(x = c2, ylim = c(0, ymax))
ggplot2::ggsave(filename = 'images/empir_p_scan_phen3.pdf', height = 3, width = 9)

plot(x = d2, ylim = c(0, ymax))
ggplot2::ggsave(filename = 'images/empir_p_scan_phen4.pdf', height = 3, width = 9)


# a1r <- scanonevar.perm(cross = my.cross,
# 											 mean.formula = 'phenotype1 ~ sex + mean.QTL.add + mean.QTL.dom',
# 											 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 											 n.perms = n.perms)
# a1p <- scanonevar.to.p.values(scanonevar = a1, perm.scan.maxes = a1r)
#
#
# b1 <- scanonevar(cross = my.cross,
# 								 mean.formula = 'phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom',
# 								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 								 return.models = TRUE)$varscan
# b1r <- scanonevar.perm(cross = my.cross,
# 											 mean.formula = 'phenotype2 ~ sex + mean.QTL.add + mean.QTL.dom',
# 											 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 											 n.perms = n.perms)
# b1p <- scanonevar.to.p.values(scanonevar = b1, perm.scan.maxes = b1r)
#
#
# c1 <- scanonevar(cross = my.cross,
# 								 mean.formula = 'phenotype3 ~ sex + mean.QTL.add + mean.QTL.dom',
# 								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 								 return.models = TRUE)$varscan
# c1r <- scanonevar.perm(cross = my.cross,
# 											 mean.formula = 'phenotype3 ~ sex + mean.QTL.add + mean.QTL.dom',
# 											 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 											 n.perms = n.perms)
# c1p <- scanonevar.to.p.values(scanonevar = c1, perm.scan.maxes = c1r)
#
#
# d1 <- scanonevar(cross = my.cross,
# 								 mean.formula = 'phenotype4 ~ sex + mean.QTL.add + mean.QTL.dom',
# 								 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 								 return.models = TRUE)$varscan
# d1r <- scanonevar.perm(cross = my.cross,
# 											 mean.formula = 'phenotype4 ~ sex + mean.QTL.add + mean.QTL.dom',
# 											 var.formula = '~sex + var.QTL.add + var.QTL.dom',
# 											 n.perms = n.perms)
# d1p <- scanonevar.to.p.values(scanonevar = d1, perm.scan.maxes = d1r)
#
#
#
# # pdf(file = '../2016_G3_PackageVQTL_Corty/images/p_scans.pdf', width = 6, height = 8)
# par(mfrow = c(4, 1), mar = c(3.1, 3.1, 3.1, 2.1))
# ylim <- c(0, 12)
# plot(x = a1p, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = ylim, title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
# plot(x = b1p, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = ylim, title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
# plot(x = c1p, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = ylim, title.cex = 1.2, line.width = 1.5, suppress.chromosome = TRUE)
# plot(x = d1p, show.equations = FALSE, legend.pos = 'right', legend.ncol = 1, ylim = ylim, title.cex = 1.2, line.width = 1.5)
# # dev.off()
