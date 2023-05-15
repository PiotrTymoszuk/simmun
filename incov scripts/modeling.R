# Robust linear modeling of metabolites as a function of infection status
# and levels of inflammatory cytokines

  insert_head()

# container ------

  incov_mod <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## responses and explanatory variables
  ## the explanatory model selection is based on
  ## the results of SIMMUN modeling by backward elimination LM
  ## and include age, infection status and inflammatory parameters

  incov_mod$response_lexicon <-
    c('tryptophan',
      'kynurenine',
      'quinolinate',
      'serotonin',
      'phenylalanine',
      'tyrosine',
      'dopamine 3-O-sulfate') %>%
    exchange(dict = globals$incov_lexicon) %>%
    compress(names_to = 'variable',
             values_to = 'label')

  incov_mod$variable_lexicon <-
    c('timepoint' = 'SARS-CoV-2',
      'age' = 'Age',
      'IL6_INF' = 'IL6',
      'IL10_INF' = 'IL10',
      'TNF_INF' = 'TNF',
      'IFNG_INF' = 'IFNG') %>%
    compress(names_to = 'variable',
             values_to = 'label')

  ## analysis table
  ## normalizing numeric variables

  incov_mod$analysis_tbl <- incov$analysis_tbl %>%
    select(all_of(incov_mod$response_lexicon$variable),
           all_of(incov_mod$variable_lexicon$variable)) %>%
    map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x) %>%
    set_names(make.names(names(.)))

  incov_mod[c("response_lexicon", "variable_lexicon")] <-
    incov_mod[c("response_lexicon", "variable_lexicon")] %>%
    map(mutate,
        variable = make.names(variable))

  ## model formulas

  incov_mod$formulas <- incov_mod$response_lexicon$variable %>%
    map(~paste(.x,
               paste(incov_mod$variable_lexicon$variable, collapse = ' + '),
               sep = ' ~ ')) %>%
    map(as.formula) %>%
    set_names(incov_mod$response_lexicon$variable)

# Construction of models -----

  insert_msg('Construction of models')

  incov_mod$models <- incov_mod$formulas %>%
    map(~make_lm(data = incov_mod$analysis_tbl,
                 formula = .x,
                 mod_fun = rlm,
                 family = NULL,
                 psi = psi.huber))

# cross-validation via caret ------

  insert_msg('Cross-validation via caret')

  set.seed(123456)

  registerDoParallel(cores = 7)

  incov_mod$caret_models <- incov_mod$formulas %>%
    map(~train(form = .x,
               data = incov_mod$analysis_tbl %>%
                 mutate(timepoint = car::recode(timepoint,
                                                "'sub-acute' = 'sub'")),
               method = 'rlm',
               metric = 'RMSE',
               trControl = trainControl(method = 'repeatedcv',
                                        number = 10,
                                        savePredictions = 'final',
                                        returnData = TRUE,
                                        returnResamp = 'final'),
               tuneGrid = data.frame(intercept = TRUE,
                                     psi = 'psi.huber'))) %>%
    map(as_caretx)

  stopImplicitCluster()

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  plan('multisession')

  incov_mod$cv_stats <- incov_mod$caret_models %>%
    extract_caret_stats

  plan('sequential')

# Plots with the R-squares and RMSE for train and CV data ------

  insert_msg('Plots of the train and CV fit stats')

  incov_mod$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(incov_mod$cv_stats,
                     subtitle = 'Robust linear modeling',
                     titles = c('Fit error, INCOV',
                                'Explained variance, INCOV')) %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(incov_mod$caret_models)),
                           labels = exchange(names(incov_mod$caret_models),
                                             dict = incov_mod$response_lexicon,
                                             key = 'variable',
                                             value = 'label')))

# Inference --------

  insert_msg('Inference and Forest plots of the estimates')

  ## inference, confidence intervals are calculated based on T distribution

  incov_mod$inference <- incov_mod$models %>%
    map(summary,
        'inference',
        ci_method = 'distribution') %>%
    map(mutate,
        var_lab = exchange(variable,
                           dict = incov_mod$variable_lexicon))

  ## Forest plots

  incov_mod$forest_plots <-
    list(x = incov_mod$inference,
         plot_title = exchange(names(incov_mod$inference),
                               dict = incov_mod$response_lexicon),
         x_lab = exchange(names(incov_mod$inference),
                          dict = incov_mod$response_lexicon)) %>%
    pmap(plot_forest,
         variable = 'var_lab',
         hide_baseline = TRUE,
         estimate_size = 2.5,
         estimate_hjust = 0.5,
         cust_theme = globals$common_theme,
         plot_subtitle = paste('n =', nrow(incov_mod$analysis_tbl)))

# Appending the Forest plots with the baselines -------

  insert_msg('Adding the baseline')

  incov_mod$intercepts <-
    incov_mod$inference[names(incov_mod$forest_plots)] %>%
    map_dfr(filter,
            level == 'baseline') %>%
    mutate(baseline_lab = paste0('baseline: ', signif(estimate, 2))) %>%
    .$baseline_lab

  incov_mod$forest_plots <-
    map2(incov_mod$forest_plots,
         incov_mod$intercepts,
         ~.x +
           labs(x = paste(.x$labels$x,
                          .y,
                          sep = ', ')))

# Caching the results -------

  insert_msg('Caching the results')

  save(incov_mod, file = './cache/incov_mod.RData')

# END -------

  insert_tail()
