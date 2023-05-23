# Supplementary Figures

  insert_head()

# container ------

  suppl <- list()

# Figure S1: consistency of psychometric tools ------

  insert_msg('Figure S1: consistency of psychometric tools')

  suppl$consistency <- cons$fa_loading_plots %>%
    map(~.x +
          labs(title = paste0(.x$labels$title, ', SIMMUN')) +
          theme(plot.tag = element_blank())) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(plot_grid(cons$omega_plot +
                          labs(title = paste0(cons$omega_plot$labels$title,
                                              ', SIMMUN')),
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

# Figure S3: univariate analysis for TRP/KYN ------

  insert_msg('Figure S3: Univariate analysis for TRP/KYN')

  ## upper panel: correlogram

  suppl$uni_trp$upper <-
    plot_grid(tst_trp$correlation$correlogram +
                theme(plot.subtitle = element_text(size = 7.4),
                      legend.position = 'none') +
                labs(subtitle = paste('Pearson correlation, n =',
                                      nrow(stigma$analysis_tbl)),
                     title = paste0(tst_trp$correlation$correlogram$labels$title,
                                    ', SIMMUN')),
              get_legend(tst_trp$correlation$correlogram),
              ncol = 2,
              rel_widths = c(0.65, 0.35),
              align = 'hv',
              axis = 'tblr')

  ## middle panel: correlations for neopterin

  suppl$uni_trp$middle <-
    tst_trp$correlation$plots[c("log_neo|trp",
                                "log_neo|log_kyn",
                                "log_neo|log_kyn_trp")] %>%
    map(~.x +
          labs(title = paste0(.x$labels$title,
                              ', SIMMUN')) +
          theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: comparisons

  suppl$uni_trp$bottom <- tst_trp$comparison$plots %>%
    map(~.x +
          labs(title = paste0(.x$labels$title, ', SIMMUN')) +
          theme(plot.subtitle = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  suppl$uni_trp <- plot_grid(suppl$uni_trp$upper,
                            suppl$uni_trp$middle,
                            suppl$uni_trp$bottom,
                            nrow = 3,
                            rel_heights = c(0.75, 1, 2),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'figure_s3_univariate_analysis_trp',
              ref_name = 'uni_trp',
              caption = paste('Effects of age, serum inflammatory markers',
                              'neopterin and neutrophil-lymphocyte ratio,',
                              'stress, depreression and SARS-CoV-2 infection',
                              'on tryptophan, kynurenine',
                              'and kynurenine/tryptophan ratio',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 230)

# Figure S4: univariate analysis for the TYR subsystem -----

  insert_msg('Figure S4: univariate analysis for TYR')

  ## upper panel: correlogram

  suppl$uni_tyr$upper <-
    plot_grid(tst_tyr$correlation$correlogram +
                theme(plot.subtitle = element_text(size = 7.4),
                      legend.position = 'none') +
                labs(subtitle = paste('Pearson correlation, n =',
                                      nrow(stigma$analysis_tbl)),
                     title = paste0(tst_tyr$correlation$correlogram$labels$title,
                                    ', SIMMUN')),
              get_legend(tst_tyr$correlation$correlogram),
              ncol = 2,
              rel_widths = c(0.55, 0.45),
              align = 'hv',
              axis = 'tblr')

  ## middle panel, correlations for age

  suppl$uni_tyr$middle <-
    tst_tyr$correlation$plots[c("age|log_phe",
                                "age|log_tyr",
                                "age|sqrt_phe_tyr")] %>%
    map(~.x +
          labs(title = paste0(.x$labels$title, ', SIMMUN')) +
          theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel

  suppl$uni_tyr$bottom <- tst_tyr$comparison$plots %>%
    map(~.x +
          labs(title = paste0(.x$labels$title, ', SIMMUN')) +
          theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  suppl$uni_tyr <- plot_grid(suppl$uni_tyr$upper,
                            suppl$uni_tyr$middle,
                            suppl$uni_tyr$bottom,
                            nrow = 3,
                            rel_heights = c(0.75, 1, 1),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'figure_s4_univariate_analysis_tyr',
              ref_name = 'uni_tyr',
              caption = paste('Effects of age, serum inflammatory marker',
                              'neopterin and SARS-CoV-2 infection',
                              'on phenylalanine, tyrosine',
                              'and phenylalanine/tyrosine ratio',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 180)

# Figure S5: INCOV, correlograms -----

  insert_msg('Figure S5: correlograms')

  suppl$corr_incov <- cyt_corr$bubbles %>%
    map(~.x + theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(cyt_corr$bubbles[[1]]),
              ncol = 2,
              rel_widths = c(0.9, 0.1)) %>%
    as_figure(label = 'figure_s5_incov_correlations',
              ref_name = 'corr_invcov',
              caption = paste('Correlation of metabolites with inflammatory',
                              'cytokine markers in the INCOV cohort.'),
              w = 180,
              h = 180)

# Saving on the disc -----

  insert_msg('Saving the supplementary figures')

  suppl %>%
    walk(pickle,
         path = './paper/supplementary figures',
         format = 'pdf',
         device = cairo_pdf)

# END -----

  insert_tail()
