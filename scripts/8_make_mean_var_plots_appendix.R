xlims <- c(-0.3, 0.5)
ylims <- c(0.6, 1.6)

plot(p1 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype1x',
																		 focal.covariate.names = 'D1M1',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p2 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype2x',
																		 focal.covariate.names = 'D1M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p3 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype3x',
																		 focal.covariate.names = 'D2M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p4 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype4',
																		 focal.covariate.names = 'D3M6',
																		 xlim = xlims,
																		 ylim = ylims))


ggplot2::ggsave(plot = p1, filename = 'images/mean_var_plot_phen1x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p2, filename = 'images/mean_var_plot_phen2x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p3, filename = 'images/mean_var_plot_phen3x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p4, filename = 'images/mean_var_plot_phen4x.pdf', height = 3, width = 4)
