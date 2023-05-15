# Main figures

  insert_head()

# container ------

  figs <- list()

# Figure 1: CONSORT -----

  insert_msg('Figure 1: study inclusion diagram')

  figs$inclusion <-
    plot_grid(ggdraw() +
                draw_image('./schemes/inclusion.png')) %>%
    as_figure(label = 'figure_1_study_inclusion_scheme',
              ref_name = 'inclusion',
              caption = paste('Inclusion scheme for the SIMMUN',
                              'and INCOV cohorts, and analysis strategy.'),
              w = 180,
              h = 180 * 2832/4044)

# Figure 2: linear modeling -----

  insert_msg('Figure 2: linear modeling')

  figs$modeling <- mod_eli$forest_plots %>%
    map(~.x + theme(legend.position = 'none'))

  figs$modeling$kyn <- figs$modeling$kyn +
    expand_limits(x = -0.15)

  figs$modeling$tyr <- figs$modeling$tyr +
    expand_limits(x = -0.13)

  figs$modeling <-
    c(figs$modeling["trp"],
      list(get_legend(figs$modeling$kyn_trp +
                        theme(legend.title = element_blank()))),
      figs$modeling[c("kyn", "kyn_trp", "tyr", "phe_tyr")]) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr',
              rel_heights = c(1, 1.35, 1),
              labels = c('A', '',
                         '', '',
                         'B', ''),
              label_size = 10) %>%
    as_figure(label = 'figure_2_linear_modeling',
              ref_name = 'modeling',
              caption = paste('Results of multi-parameter linear modeling',
                              'of tryptophan,',
                              'kynurenine, kynurenine/tryptophan ratio,',
                              'tyrosine and phenylalanine/tyrosine ratio',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 180)

# Figure 3: univariate analysis for TRP/KYN ------

  insert_msg('Figure 3: Univariate analysis for TRP/KYN')

  ## upper panel: correlogram

  figs$uni_trp$upper <-
    plot_grid(tst_trp$correlation$correlogram +
                theme(plot.subtitle = element_text(size = 7.4),
                      legend.position = 'none') +
                labs(subtitle = paste('Pearson correlation, n =',
                                      nrow(stigma$analysis_tbl))),
              get_legend(tst_trp$correlation$correlogram),
              ncol = 2,
              rel_widths = c(0.65, 0.35),
              align = 'hv',
              axis = 'tblr')

  ## middle panel: correlations for neopterin

  figs$uni_trp$middle <-
    tst_trp$correlation$plots[c("log_neo|trp",
                                "log_neo|log_kyn",
                                "log_neo|log_kyn_trp")] %>%
    map(~.x + theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel: comparisons

  figs$uni_trp$bottom <- tst_trp$comparison$plots %>%
    map(~.x + theme(plot.subtitle = element_text(size = 7))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  figs$uni_trp <- plot_grid(figs$uni_trp$upper,
                            figs$uni_trp$middle,
                            figs$uni_trp$bottom,
                            nrow = 3,
                            rel_heights = c(0.75, 1, 2),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'figure_3_univariate_analysis_trp',
              ref_name = 'uni_trp',
              caption = paste('Effects of age, serum inflammatory markers',
                              'neopterin and neutrophil-lymphocyte ratio,',
                              'stress, depreression and SARS-CoV-2 infection',
                              'on tryptophan, kynurenine',
                              'and kynurenine/tryptophan ratio',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 230)

# Figure 4: univariate analysis for the TYR subsystem -----

  insert_msg('Figure 4: univariate analysis for TYR')

  ## upper panel: correlogram

  figs$uni_tyr$upper <-
    plot_grid(tst_tyr$correlation$correlogram +
                theme(plot.subtitle = element_text(size = 7.4),
                      legend.position = 'none') +
                labs(subtitle = paste('Pearson correlation, n =',
                                      nrow(stigma$analysis_tbl))),
              get_legend(tst_tyr$correlation$correlogram),
              ncol = 2,
              rel_widths = c(0.55, 0.45),
              align = 'hv',
              axis = 'tblr')

  ## middle panel, correlations for age

  figs$uni_tyr$middle <-
    tst_tyr$correlation$plots[c("age|log_phe",
                                "age|log_tyr",
                                "age|sqrt_phe_tyr")] %>%
    map(~.x +
          theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  ## bottom panel

  figs$uni_tyr$bottom <- tst_tyr$comparison$plots %>%
    map(~.x + theme(plot.subtitle = element_text(size = 7.4))) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr')

  figs$uni_tyr <- plot_grid(figs$uni_tyr$upper,
                            figs$uni_tyr$middle,
                            figs$uni_tyr$bottom,
                            nrow = 3,
                            rel_heights = c(0.75, 1, 1),
                            labels = LETTERS,
                            label_size = 10) %>%
    as_figure(label = 'figure_4_univariate_analysis_tyr',
              ref_name = 'uni_tyr',
              caption = paste('Effects of age, serum inflammatory marker',
                              'neopterin and SARS-CoV-2 infection',
                              'on phenylalanine, tyrosine',
                              'and phenylalanine/tyrosine ratio',
                              'in the SIMMUN cohort.'),
              w = 180,
              h = 180)

# Figure 5: multi-parameter modeling of neurotransmitters in the INCOV cohort -----

  insert_msg('Figure 5: INCOV neurotransmitters')

  ## upper panel: fit stats

  figs$incov_neuro$upper <- incov_neuro$fit_stat_plots %>%
    map(~.x +
          labs(subtitle = 'Robust linear modeling') +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(incov_neuro$fit_stat_plots[[1]] +
                           theme(legend.position = 'bottom')),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## bottom panel: Forest plots with inference measures
  ## hacking for re-ordering the Y axis

  figs$incov_neuro$plot_order <-
    c('age' = 'age',
      'timepointacute' = 'SARS-CoV-2',
      'timepointsub-acute' = 'SARS-CoV-2',
      'timepointrecovery' = 'SARS-CoV-2',
      'tryptophan' = 'metabolites',
      'kynurenine' = 'metabolites',
      'quinolinate' = 'metabolites',
      'phenylalanine' = 'metabolites',
      'tyrosine' = 'metabolites',
      'IL6_INF' = 'inflammation',
      'IL10_INF' = 'inflammation',
      'TNF_INF' = 'inflammation',
      'IFNG_INF' = 'inflammation') %>%
    compress(names_to = 'parameter',
             values_to = 'subset') %>%
    mutate(plot_order = 1:nrow(.),
           subset = factor(subset,
                           rev(c('inflammation',
                                 'metabolites',
                                 'SARS-CoV-2',
                                 'age'))))

  figs$incov_neuro$bottom <- incov_neuro$forest_plots

  for(i in names(figs$incov_neuro$bottom)) {

    figs$incov_neuro$bottom[[i]]$data <-
      figs$incov_neuro$bottom[[i]]$data %>%
      select(-plot_order) %>%
      left_join(figs$incov_neuro$plot_order,
                by = 'parameter') %>%
      mutate(y_ax = stri_replace(y_ax,
                                 fixed = 'SARS-CoV-2: ',
                                 replacement = ''))

  }

  figs$incov_neuro$bottom <- figs$incov_neuro$bottom %>%
    map(~.x +
          expand_limits(x = -0.5) +
          expand_limits(x = 1.11) +
          facet_grid(subset ~ .,
                     scales = 'free',
                     space = 'free') +
          theme(legend.position = 'none')) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              align = 'hv',
              axis = 'tblr') %>%
    plot_grid(get_legend(figs$incov_neuro$bottom[[1]] +
                           theme(legend.position = 'bottom',
                                 legend.title = element_blank())),
              nrow = 2,
              rel_heights = c(0.9, 0.1))

  ## the entire figure

  figs$incov_neuro <-
    plot_grid(figs$incov_neuro$upper,
              figs$incov_neuro$bottom,
              nrow = 2,
              rel_heights = c(0.3, 0.7),
              labels = LETTERS,
              label_size = 10) %>%
    as_figure(label = 'figure_5_incov_neurotransmitter_modeling',
              ref_name = 'incov_neuro',
              caption = paste('Results of multi-parameter robust',
                              'linear modeling of serum levels of',
                              'serotinin and dopamine sulfate',
                              'in the INCOV cohort.'),
              w = 180,
              h = 210)

# Figure 6: INCOV, time course ----

  insert_msg('Figure 6: INCOV time course')

  figs$time_course <- names(time_rlm$ribbon_plots) %>%
    map(~list(time_rlm$ribbon_plots[[.x]] +
                theme(legend.position = 'none'),
              time_rlm$forest_plots[[.x]] +
                theme(legend.position = 'none') +
                labs(title = ' '),
              get_legend(time_rlm$forest_plots[[.x]] +
                           theme(legend.title = element_blank())))) %>%
    reduce(c) %>%
    plot_grid(plotlist = .,
              ncol = 3,
              align = 'hv',
              axis = 'tblr',
              rel_widths = c(1, 1, 0.3),
              labels = c('A', '', '',
                         'B', '', '',
                         'C', '', ''),
              label_size = 10) %>%
    as_figure(label = 'figure_6_incov_time_course',
              ref_name = 'time_ourse',
              caption = paste('Time course of inflammatory cytokines',
                              'tryptophan, tyrosine and their metabolites',
                              'during SARS-CoV infection and recovery',
                              'in the INCOV cohort.'),
              w = 180,
              h = 210)

# Figure 7: INCOV, correlation -----

  insert_msg('Figure 7: INCOV correlations')

  figs$incov_corr <- incov_mds$network_plots %>%
    map(~.x +
          theme(legend.position = 'none',
                plot.margin = globals$common_margin)) %>%
    plot_grid(plotlist = .,
              ncol = 2,
              labels = c('A', '', '', '',
                         'B', '', '', ''),
              label_size = 10) %>%
    plot_grid(get_legend(incov_mds$network_plots$`serotonin.sub-acute`),
              ncol = 2,
              rel_widths = c(0.88, 0.12)) %>%
    as_figure(label = 'figure_7_incov_correlation',
              ref_name = 'incov_corr',
              caption = paste('Correlation of serum levels of inflammatory',
                              'cytokines, serotonin and domamine sulfate,',
                              'their precursors and competitor pathway products',
                              'in healthy controls and SARS-CoV-2',
                              'patients of the INCOV cohort.'),
              w = 180,
              h = 220)

# Saving the figures on the disc ------

  insert_msg('Saving figures on the disc')

  figs %>%
    walk(pickle,
         format = 'pdf',
         path = './paper/figures',
         device = cairo_pdf)

# END -----

  insert_tail()
