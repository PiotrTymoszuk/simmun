# cyt_correlation of inflammatory and metabolic INCOV variables (Spearman)

  insert_head()

# container -----

  cyt_corr <- list()

# globals ------

  insert_msg('Globals')

  ## analysis table

  cyt_corr$analysis_tbl <- incov$analysis_tbl %>%
    select( patient_id,
            timepoint,
            any_of(globals$incov_lexicon$variable)) %>%
    blast(timepoint)

  ## variable pairs

  cyt_corr$pairs <- globals$incov_proteins %>%
    map(function(prot) globals$incov_metabolites %>%
          map(~c(prot, .x))) %>%
    unlist(recursive = FALSE)

# Correlation analysis ------

  insert_msg('Correlation')

  cyt_corr$test <- cyt_corr$analysis_tbl %>%
    map(function(dat) cyt_corr$pairs %>%
          map_dfr(~correlate_variables(dat,
                                       variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       ci = TRUE)))
  ## FDR cyt_correction, plot elements

  cyt_corr$test <- cyt_corr$test %>%
    compress(names_to = 'timepoint') %>%
    re_adjust %>%
    mutate(timepoint = factor(timepoint, names(cyt_corr$test)),
           plot_lab = paste0(signif(estimate, 2),
                             ' [', signif(lower_ci, 2),
                             ' - ', signif(upper_ci, 2), ']'),
           plot_cap = paste0(plot_lab, ', ', significance, ', n = ', n),
           x_lab = exchange(variable1,
                            dict = globals$incov_lexicon),
           x_lab = paste0(x_lab, ', rank'),
           y_lab = exchange(variable2,
                            dict = globals$incov_lexicon),
           y_lab = paste0(y_lab, ', rank')) %>%
    blast(timepoint)

# Bubble plots -------

  insert_msg('Bubble plots')

  cyt_corr$bubbles <-
    list(data = cyt_corr$test %>%
           map(mutate,
               variable2 = factor(variable2,
                                  c('phenylalanine',
                                    'tyrosine',
                                    'dopamine 3-O-sulfate',
                                    'quinolinate',
                                    'kynurenine',
                                    'tryptophan',
                                    'serotonin')),
               variable1 = factor(variable1,
                                  c('IL6_INF',
                                    'IL10_INF',
                                    'TNF_INF',
                                    'IFNG_INF'))),
         plot_title =  c('uninfected, INCOV',
                         'acute SARS-CoV-2 infection, INCOV',
                         'sub-acute SARS-CoV-2 infection, INCOV',
                         'SARS-CoV-2 infection recovery, INCOV')) %>%
    pmap(corr_buble,
         signif_only = FALSE,
         rotate_x_labs = FALSE)


  ## faceting for easier interpretation

  for(i in names(cyt_corr$bubbles)) {

    cyt_corr$bubbles[[i]]$data <-
      cyt_corr$bubbles[[i]]$data %>%
      mutate(pathway = ifelse(variable2 %in% c('tryptophan',
                                               'kynurenine',
                                               'quinolinate',
                                               'serotonin'),
                              'serotonin/kynurenine', 'catecholamine'),
             pathway = factor(pathway, c('serotonin/kynurenine', 'catecholamine')))

  }

  cyt_corr$bubbles <- cyt_corr$bubbles %>%
    map(~.x +
          facet_grid(pathway ~ ., space = 'free', scales = 'free') +
          scale_y_discrete(labels = function(x) exchange(x, dict = globals$incov_lexicon)) +
          scale_x_discrete(labels = function(x) exchange(x, dict = globals$incov_lexicon)))
# Classical correlation plots -------

  insert_msg('Classical correlation plots')


  ## plotting the ranks, as suggested by a Reviewer

  cyt_corr$point_plots <- list(dat = cyt_corr$analysis_tbl %>%
                             map(~map_dfc(.x, rank)),
                           stat = cyt_corr$test,
                           tit = paste('INCOV,', names(cyt_corr$analysis_tbl))) %>%
    pmap(function(dat, stat, tit) list(variables = cyt_corr$pairs,
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
           set_names(map(cyt_corr$pairs, paste, collapse = '_')))

# END ------

  insert_tail()
