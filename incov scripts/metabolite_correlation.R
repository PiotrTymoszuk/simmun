# Inter-correlation of neurotransmitter-related metabolites
# In the INCOV cohort

  insert_head()

# container ------

  met_corr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variables

  met_corr$variables <-
    c('tryptophan',
      'kynurenine',
      'quinolinate',
      'serotonin',
      'phenylalanine',
      'tyrosine',
      'dopamine 3-O-sulfate')

  met_corr$lexicon <- globals$incov_lexicon %>%
    filter(variable %in% met_corr$variables) %>%
    mutate(variable = make.names(variable))

  ## analysis tables

  met_corr$analysis_tbl <- incov$analysis_tbl %>%
    select(timepoint,
           all_of(met_corr$variables)) %>%
    set_names(make.names(names(.))) %>%
    blast(timepoint) %>%
    map(select, -timepoint)

  met_corr$variables <- make.names(met_corr$variables)

  ## variable pairs

  met_corr$pairs <- met_corr$variables %>%
    combn(m = 2, simplify = FALSE)

# Correlation analysis: Spearman's rank test -----

  insert_msg('Correlation')

  met_corr$test <- met_corr$analysis_tbl %>%
    map(function(data) met_corr$pairs %>%
          map_dfr(~correlate_variables(data,
                                       variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       ci = TRUE,
                                       pub_styled = FALSE))) %>%
    compress(names_to = 'timepoint') %>%
    re_adjust %>%
    mutate(timepoint = factor(timepoint, names(met_corr$analysis_tbl)),
           est_lab = paste0(signif(estimate, 2), ' [',
                            signif(lower_ci, 2), ' - ',
                            signif(upper_ci, 2), ']'),
           plot_cap = paste(est_lab, significance, sep = ', '),
           x_lab = exchange(variable1,
                            dict = globals$incov_lexicon),
           x_lab = paste0(x_lab, ', rank'),
           y_lab = exchange(variable2,
                            dict = globals$incov_lexicon),
           y_lab = paste0(y_lab, ', rank')) %>%
    blast(timepoint)

# Correlograms ------

  insert_msg('Correlograms')

  met_corr$bubbles <-
    list(data = met_corr$test,
         plot_title =  c('uninfected, INCOV',
                         'acute SARS-CoV-2 infection, INCOV',
                         'sub-acute SARS-CoV-2 infection, INCOV',
                         'SARS-CoV-2 infection recovery, INCOV')) %>%
    pmap(corr_buble,
         signif_only = FALSE,
         rotate_x_labs = FALSE)

  ## extra adjustment of scales

  met_corr$bubbles <- met_corr$bubbles %>%
    map(~.x +
          scale_x_discrete(limits = met_corr$variables,
                           labels = function(x) exchange(x, dict = met_corr$lexicon)) +
          scale_y_discrete(limits = met_corr$variables,
                           labels = function(x) exchange(x, dict = met_corr$lexicon)))

# Classical correlation plots -------

  insert_msg('Classical correlation plots')

  ## visualized are ranks, as requested by a Reviewer

  met_corr$point_plots <- list(dat = met_corr$analysis_tbl %>%
                                 map(~map_dfc(.x, rank)),
                               stat = met_corr$test,
                               tit = paste('INCOV,', names(met_corr$analysis_tbl))) %>%
    pmap(function(dat, stat, tit) list(variables = met_corr$pairs,
                                       plot_subtitle = stat$plot_cap,
                                       x_lab = stat$x_lab,
                                       y_lab = stat$y_lab) %>%
           pmap(plot_correlation,
                data = dat,
                type = 'correlation',
                point_hjitter = 0,
                point_wjitter = 0,
                point_color = 'coral3',
                plot_title = tit,
                show_trend = FALSE,
                cust_theme = globals$common_theme) %>%
           map(~.x +
                 geom_smooth(method = 'lm') +
                 theme(plot.tag = element_blank())) %>%
           set_names(map(met_corr$pairs, paste, collapse = '_')))

# END ------

  insert_tail()
