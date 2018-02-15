#### SIMULATE CROSS AND PHENOTYPES ####
library(tidyr); library(dplyr); library(qtl); library(vqtl);

r <- readRDS('scans_2017-05-07 12:09:54.RDS')

for (phen.name in names(r[['sos']])) {

	p <- plot(x = r[['sovs']][[phen.name]], y = r[['sos']][[phen.name]], plotting.units = 'LOD') +
		ggplot2::coord_cartesian(ylim = c(0, 5))

	ggplot2::ggsave(plot = p, height = 2.5, width = 9,
									filename = paste0('images/LOD_scan_', phen.name,'.pdf'))

	p <- plot(x = r[['sovs']][[phen.name]], y = r[['sos']][[phen.name]]) +
		ggplot2::coord_cartesian(ylim = c(0, 2.5))

	ggplot2::ggsave(plot = p, height = 2.5, width = 9,
									filename = paste0('images/empir_p_scan_', phen.name,'.pdf'))
}


for (phen.name in names(r[['soxs']])) {

	p <- plot(x = r[['sovxs']][[phen.name]], y = r[['soxs']][[phen.name]], plotting.units = 'LOD') +
		ggplot2::coord_cartesian(ylim = c(0, 6))

	ggplot2::ggsave(plot = p, height = 2.5, width = 9,
									filename = paste0('images/LOD_scan_', phen.name,'.pdf'))

	p <- plot(x = r[['sovxs']][[phen.name]], y = r[['soxs']][[phen.name]]) +
		ggplot2::coord_cartesian(ylim = c(0, 2.5))

	ggplot2::ggsave(plot = p, height = 2.5, width = 9,
									filename = paste0('images/empir_p_scan_', phen.name,'.pdf'))
}
