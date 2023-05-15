# Supplementary Figures

  insert_head()

# container ------

  suppl <- list()

# Figure S1: consistency of psychometric tools ------

  insert_msg('Figure S1: consistency of psychometric tools')

  suppl$consistency <- cons$fa_loading_plots %>%
    map(~.x + theme(plot.tag = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(plot_grid(cons$omega_plot,
                        ncol = 2),
              nrow = 2,
              rel_heights = c(1, 0.7),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'figure_s1_psycho_consistency',
              ref_name = 'consistency',
              caption = paste('Consistency of the PSS4 stress',
                              'HADS depression and HADS anxiety psychometric',
                              'tools in the SIMMUN cohort.'),
              w = 180,
              h = 140)

# Figure S2: error and R-squared, multiple linear regression ------

  insert_msg('Figure S2: Fits stats, multiple regression')

  suppl$fit_stats <- mod_eli$fit_stat_plots %>%
    map2(c('pseudo_log', 'identity'),
         ~.x +
           scale_x_continuous(trans = .y) +
           theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              labels = LETTERS,
              label_size = 10) %>%
    plot_grid(get_legend(mod_eli$fit_stat_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_s2_fit_stats',
              ref_name = 'fit_stats',
              caption = paste('Root mean square error and R-squared',
                              'statistics for multi-parameter linear models',
                              'of trayptophan, tyrosine and their metabolites',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 110)

# Figure S3: scatter plots for TRP, KYN and KYN/TRP -----

  insert_msg('Figure S3: scatter plots for TRP, KYN and KYN/TRP')

  suppl$corr_trp <-
    tst_trp$correlation$plots[c("age|trp", "age|log_kyn", "age|log_kyn_trp",
                                "log_nlr|trp", "log_nlr|log_kyn",
                                "log_nlr|log_kyn_trp", "pss_stress_score|trp",
                                "pss_stress_score|log_kyn",
                                "pss_stress_score|log_kyn_trp")] %>%
    map(~.x + theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr',
              labels = c('A', '', '',
                         'B', '', '',
                         'C', '', ''),
              label_size = 10) %>%
    as_figure(label = 'figure_s3_trp_system_correlation',
              ref_name = 'corr_trp',
              caption = paste('Correlation of tryptophan and',
                              'its metabolites with age, serum levels',
                              'of neopterin, neutrophil - lymphocyte ratio,',
                              'and mental stress scoring in the SIMMUN',
                              'cohort.'),
              w = 180,
              h = 180)

# Figure S4: scatter plots for PHE and TYR -----

  insert_msg('Figure S4: scatter plots for phe and tyr')

  suppl$corr_tyr <-
    tst_tyr$correlation$plots[c("log_neo|log_phe",
                                "log_neo|log_tyr",
                                "log_neo|sqrt_phe_tyr")] %>%
    map(~.x + theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr') %>%
    as_figure(label = 'figure_S4_tyr_system_correlation',
              ref_name = 'corr_tyr',
              caption = paste('Correlation of phenylalanine',
                              'and tyrosine with serum levels',
                              'of neopterin in the SIMMUN cohort.'),
              w = 180,
              h = 70)

# Saving on the disc -----

  insert_msg('Saving the supplementary figures')

  suppl %>%
    walk(pickle,
         path = './paper/supplementary figures',
         format = 'pdf',
         device = cairo_pdf)

# END -----

  insert_tail()
