#### MEAN_VAR PLOTS FOR BODY OF PAPER ####
xlims <- c(-0.5, 0.3)
ylims <- c(0.6, 1.5)

plot(p1 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype1',
																		 focal.groups = 'D1M1',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p2 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype2',
																		 focal.groups = 'D1M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p3 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype3',
																		 focal.groups = 'D2M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p4 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype4',
																		 focal.groups = 'D3M6',
																		 xlim = xlims,
																		 ylim = ylims))


ggplot2::ggsave(plot = p1, filename = 'images/mean_var_plot_phen1.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p2, filename = 'images/mean_var_plot_phen2.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p3, filename = 'images/mean_var_plot_phen3.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p4, filename = 'images/mean_var_plot_phen4.pdf', height = 3, width = 4)


#### MEAN VAR PLOTS FOR APPENDIX ####
xlims <- c(-0.4, 0.5)
ylims <- c(0.6, 1.6)

plot(p1 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype1x',
																		 focal.groups = 'D1M1',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p2 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype2x',
																		 focal.groups = 'D1M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p3 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype3x',
																		 focal.groups = 'D2M6',
																		 xlim = xlims,
																		 ylim = ylims))

plot(p4 <- mean_var_plot_model_based(cross = my.cross,
																		 phenotype.name = 'phenotype4',
																		 focal.groups = 'D3M6',
																		 xlim = xlims,
																		 ylim = ylims))


ggplot2::ggsave(plot = p1, filename = 'images/mean_var_plot_phen1x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p2, filename = 'images/mean_var_plot_phen2x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p3, filename = 'images/mean_var_plot_phen3x.pdf', height = 3, width = 4)
ggplot2::ggsave(plot = p4, filename = 'images/mean_var_plot_phen4x.pdf', height = 3, width = 4)
