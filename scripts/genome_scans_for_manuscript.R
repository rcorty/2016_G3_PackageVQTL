#### SIMULATE CROSS AND PHENOTYPES ####
library(tidyr); library(dplyr); library(qtl); library(vqtl);

set.seed(27599)

# simulate the cross
my.cross <- qtl::sim.cross(map = qtl::sim.map(len = rep(100, 3), n.mar = 11, eq.spacing = TRUE, include.x = FALSE, anchor.tel = TRUE),
													 n.ind = 400,
													 type = 'f2')
my.cross$pheno$batch <- sample(x = 10, size = 400, replace = TRUE)
my.cross <- qtl::calc.genoprob(cross = my.cross, step = 2)
batch.effects <- runif(n = 10, min = -0.5, max = 0.5)

# simulate phenotype1 through phenotype4, for the body of the paper
my.cross$pheno$phenotype1 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0,
																	 sd = 1)
my.cross$pheno$phenotype2 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0.31*(my.cross$geno$`1`$data[,6] - 2),
																	 sd = 1)
my.cross$pheno$phenotype3 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0,
																	 sd = exp(0.2*(my.cross$geno$`2`$data[,6] - 2)))
my.cross$pheno$phenotype4 <- rnorm(n = qtl::nind(my.cross),
																	 mean = 0.23*(my.cross$geno$`3`$data[,6] - 2),
																	 sd = exp(0.17*(my.cross$geno$`3`$data[,6] - 2)))

# simulate phenotype1x through phenotype4x, for the appendix
my.cross$pheno$phenotype1x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0,
																		sd = exp(0.5*batch.effects[my.cross$pheno$batch]))
my.cross$pheno$phenotype2x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0.45*(my.cross$geno$`1`$data[,6] - 2),
																		sd = exp(batch.effects[my.cross$pheno$batch]))
my.cross$pheno$phenotype3x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0,
																		sd = exp(0.24*(my.cross$geno$`2`$data[,6] - 2) + (batch.effects[my.cross$pheno$batch])))
my.cross$pheno$phenotype4x <- rnorm(n = qtl::nind(my.cross),
																		mean = 0.22*(my.cross$geno$`3`$data[,6] - 2),
																		sd = exp(0.15*(my.cross$geno$`3`$data[,6] - 2) + (batch.effects[my.cross$pheno$batch])))


sos <- sovs <- list()
for (focal.phen.name in c('phenotype1', 'phenotype2', 'phenotype3', 'phenotype4')) {
	##### LOD SCORE GENOME SCAN ####
	so <- qtl::scanone(cross = my.cross,
										 pheno.col = focal.phen.name,
										 addcovar = my.cross$pheno$batch)
	sov <- scanonevar(cross = my.cross,
										mean.formula = formula(paste(focal.phen.name, '~mean.QTL.add + mean.QTL.dom')),
										var.formula = ~ var.QTL.add + var.QTL.dom)

	#### PERMUTATIONS FOR EMPIRICAL P-VALUES ####
	perms <-  qtl::scanone(cross = my.cross,
												 pheno.col = focal.phen.name,
												 n.perm = 1000,
												 verbose = FALSE)
	the.evd <- evd::fgev(x = perms)
	so$empir.p <- evd::pgev(q = so$lod,
													loc = fitted(the.evd)[1],
													scale = fitted(the.evd)[2],
													shape = fitted(the.evd)[3],
													lower.tail = FALSE)

	sov <- scanonevar.perm(sov = sov, n.perms = 1000, random.seed = 27599, n.cores = 40)

	sos[[focal.phen.name]] <- so
	sovs[[focal.phen.name]] <- sov

	message('Done with permutations for ', focal.phen.name)

}



soxs <- sovxs <- list()

for (focal.phen.name in c('phenotype1x', 'phenotype2x', 'phenotype3x', 'phenotype4x')) {

	##### LOD SCORE GENOME SCAN ####
	so <- qtl::scanone(cross = my.cross,
										 pheno.col = focal.phen.name,
										 addcovar = my.cross$pheno$batch)
	sov <- scanonevar(cross = my.cross,
										mean.formula = formula(paste(focal.phen.name, '~ batch + mean.QTL.add + mean.QTL.dom')),
										var.formula = ~batch + var.QTL.add + var.QTL.dom)
	message('Done with ', focal.phen.name)


	#### PERMUTATIONS FOR EMPIRICAL P-VALUES ####
	perms <-  qtl::scanone(cross = my.cross,
												 pheno.col = focal.phen.name,
												 addcovar = my.cross$pheno$batch,
												 n.perm = 1000,
												 verbose = FALSE)
	the.evd <- evd::fgev(x = perms)
	so$empir.p <- evd::pgev(q = so$lod,
													loc = fitted(the.evd)[1],
													scale = fitted(the.evd)[2],
													shape = fitted(the.evd)[3],
													lower.tail = FALSE)

	sov <- scanonevar.perm(sov = sov, n.perms = 1000, random.seed = 27599, n.cores = 40)

	soxs[[focal.phen.name]] <- so
	sovxs[[focal.phen.name]] <- sov

	message('Done with permutations for ', focal.phen.name)

}

results <- list(sos = sos,
								sovs = sovs,
								soxs = soxs,
								sovxs = sovxs)

saveRDS(object = results,
				file = paste0('scans_', Sys.time(),'.RDS'))
