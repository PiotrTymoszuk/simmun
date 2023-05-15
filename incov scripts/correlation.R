# Correlation analysis results to be presented in the paper
# Levels of 5-HT and DA sulfate are correlated ith normalized serum
# concentrations of the relevant metabolites and inflammatory cytokines

  insert_head()

# container -----

  corr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variables

  corr$responses <- c(serotonin = 'serotonin',
                      dopamine = 'dopamine 3-O-sulfate')

  corr$variables <-
    list(serotonin = c('tryptophan',
                       'kynurenine',
                       'quinolinate',
                       'IL6_INF',
                       'IL10_INF',
                       'TNF_INF',
                       'IFNG_INF'),
         dopamine = c('phenylalanine',
                      'tyrosine',
                      'IL6_INF',
                      'IL10_INF',
                      'TNF_INF',
                      'IFNG_INF'))

  ## analysis tables

  corr$analysis_tbl <- incov$analysis_tbl %>%
    select(timepoint,
           all_of(unname(corr$responses)),
           all_of(reduce(corr$variables, union))) %>%
    set_names(make.names(names(.))) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"),
           timepoint = factor(timepoint,
                              c('uninfected', 'acute',
                                'sub-acute', 'recovery'))) %>%
    blast(timepoint)

  ## lexicon

  corr$lexicon <- globals$incov_lexicon %>%
    mutate(variable = make.names(variable),
           class = ifelse(stri_detect(variable, regex = '_INF$'),
                          'inflammation',
                          'metabolites'))

  corr$responses <- make.names(corr$responses) %>%
    set_names(names(corr$responses))

  corr$variables <- corr$variables %>%
    map(make.names)

  ## variable pairs

  corr$pairs <-
    map2(corr$responses,
         corr$variables,
         function(x, y) map(y, ~c(x, .x))) %>%
    unlist(recursive = FALSE)

  corr$pairs <- corr$pairs %>%
    set_names(map(corr$pairs, paste, collapse = '|'))

  ## n numbers per timepoint

  corr$n_labs <- corr$analysis_tbl %>%
    map_dbl(nrow) %>%
    map2_chr(names(.), .,
             paste, sep = '\nn = ') %>%
    set_names(names(corr$analysis_tbl))

# Correlation -----

  insert_msg('Correlation')

  corr$test <- corr$analysis_tbl %>%
    map(function(data) corr$pairs %>%
          map_dfr(~correlate_variables(data, variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       ci = TRUE,
                                       pub_styled = FALSE))) %>%
    compress(names_to = 'timepoint') %>%
    re_adjust %>%
    mutate(corr_id = paste(variable1, variable2, sep = '|'),
           est_lab = paste(signif(estimate, 2), ' [',
                           signif(lower_ci, 2), ' - ',
                           signif(upper_ci, 2), ']'),
           plot_cap = paste0(est_lab, ', ', significance, ', n = ', n),
           variable1 = factor(variable1, unname(corr$responses)),
           timepoint = factor(timepoint, names(corr$analysis_tbl)),
           x_lab = exchange(variable1,
                            dict = corr$lexicon),
           x_lab = paste0(x_lab, ', rank'),
           y_lab = exchange(variable2,
                            dict = corr$lexicon),
           y_lab = paste0(y_lab, ', rank')) %>%
    left_join(corr$lexicon[c('variable', 'class')] %>%
                set_names(c('variable2', 'class')),
              by = 'variable2') %>%
    blast(variable1) %>%
    map2(., corr$variables,
         ~mutate(.x, variable2 = factor(variable2, .y))) %>%
    set_names(names(corr$responses))

# Bubble plots -------

  insert_msg('Bubble plots')

  cyt_corr$bubbles <-
    list(data = corr$test %>%
           map(mutate,
               variable1 = variable2,
               variable2 = timepoint),
         plot_title = c('5-HT, INCOV', 'DA sulfate, INCOV')) %>%
    pmap(corr_buble,
         signif_only = FALSE) %>%
    map(~.x +
          labs(subtitle = "Spearman's correlation") +
          scale_x_discrete(labels = function(x) exchange(x, dict = corr$lexicon)) +
          scale_y_discrete(labels = corr$n_labs,
                           limits = rev(names(corr$analysis_tbl))) +
          facet_grid(. ~ class,
                     space = 'free',
                     scales = 'free'))

# Classical correlation plots -------

  insert_msg('Classical correlation plots')

  ## plotting the ranks, as suggested by a Reviewer

  corr$point_plots <- list(dat = corr$analysis_tbl %>%
                             map(~map_dfc(.x, rank)),
                           stat = corr$test %>%
                             compress(names_to = 'variable1') %>%
                             mutate(timepoint = factor(timepoint,
                                                       names(corr$analysis_tbl))) %>%
                             blast(timepoint),
                           tit = paste('INCOV,', names(cyt_corr$analysis_tbl))) %>%
    pmap(function(dat, stat, tit) list(variables = corr$pairs,
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
           set_names(names(corr$pairs)))

# END ------

  insert_tail()
