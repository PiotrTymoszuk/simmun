# Modeling of neurotransmitter levels as a function of treir precursors,
# competitive decay products, infection status and inflammatory cytokines
#
# A short note on cross-validation: the model are validated by 10-fold CV
# where the folds are constructed in a timepoint-balanced wise,
# or more precisely speaking timepoint-stratified CV,
# see: tools/tools.R for code details and rationale

  insert_head()

# container ------

  incov_neuro <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## responses and explanatory variables
  ## the explanatory model selection is based on
  ## the results of SIMMUN modeling by backward elimination LM
  ## and include age, infection status and inflammatory parameters

  incov_neuro$response_lexicon <-
    c('serotonin',
      'dopamine 3-O-sulfate') %>%
    exchange(dict = globals$incov_lexicon) %>%
    compress(names_to = 'variable',
             values_to = 'label')

  incov_neuro$variable_lexicons <-
    list(serotonin = c('tryptophan' = 'TRP',
                       'kynurenine' = 'KYN',
                       'quinolinate' = 'QUIN'),
         dopamine = c('phenylalanine' = 'PHE',
                      'tyrosine' = 'TYR')) %>%
    map(~c(.x,
           c('timepoint' = 'SARS-CoV-2',
             'age' = 'Age',
             'sex' = 'Sex',
             'bmi_class' = 'BMI',
             'IL6_INF' = 'IL6',
             'IL10_INF' = 'IL10',
             'TNF_INF' = 'TNF',
             'IFNG_INF' = 'IFNG'))) %>%
    map(compress,
        names_to = 'variable',
        values_to = 'label')

  ## analysis table, expressing age in decades
  ## normalizing numeric variables

  incov_neuro$analysis_tbl <- incov$analysis_tbl %>%
    select(all_of(incov_neuro$response_lexicon$variable),
           all_of(union(incov_neuro$variable_lexicons$serotonin$variable,
                        incov_neuro$variable_lexicons$dopamine$variable))) %>%
    mutate(age = age/10) %>%
    map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x) %>%
    set_names(make.names(names(.)))

  incov_neuro$response_lexicon <- incov_neuro$response_lexicon %>%
    mutate(variable = make.names(variable))

  ## model formulas

  incov_neuro$formulas <-
    map2(incov_neuro$response_lexicon$variable,
         incov_neuro$variable_lexicons,
         ~paste(.x,
                paste(.y[['variable']], collapse = ' + '),
                sep = ' ~ ')) %>%
    map(as.formula) %>%
    set_names(names(incov_neuro$variable_lexicons))

  ## trainControl object

  set.seed(123456)

  incov_neuro$trainControl <-
    time_balanced_folds(data = incov_neuro$analysis_tbl,
                        time_variable = 'timepoint',
                        number = 10)

# Construction of models -----

  insert_msg('Construction of models')

  incov_neuro$models <- incov_neuro$formulas %>%
    map(~make_lm(data = incov_neuro$analysis_tbl,
                 formula = .x,
                 mod_fun = rlm,
                 family = NULL,
                 psi = psi.huber,
                 method = 'MM'))

# cross-validation via caret ------

  insert_msg('Cross-validation via caret')

  set.seed(123456)

  registerDoParallel(cores = 7)

  incov_neuro$caret_models <- incov_neuro$formulas %>%
    map(~train(form = .x,
               data = incov_neuro$analysis_tbl %>%
                 mutate(timepoint = car::recode(timepoint,
                                                "'sub-acute' = 'sub'")),
               method = mm_rlm,
               metric = 'RMSE',
               trControl = incov_neuro$trainControl,
               tuneGrid = data.frame(intercept = TRUE,
                                     psi = 'psi.huber'))) %>%
    map(as_caretx)

  stopImplicitCluster()

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  plan('multisession')

  incov_neuro$cv_stats <- incov_neuro$caret_models %>%
    extract_caret_stats

  plan('sequential')

# Plots with the R-squares and RMSE for train and CV data ------

  insert_msg('Plots of the train and CV fit stats')

  incov_neuro$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(incov_neuro$cv_stats,
                     subtitle = 'Robust linear modeling',
                     titles = c('Fit error, INCOV',
                                'Explained variance, INCOV')) %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(incov_neuro$caret_models)),
                           labels = exchange(names(incov_neuro$caret_models),
                                             dict = incov_neuro$response_lexicon,
                                             key = 'variable',
                                             value = 'label')))

# Inference --------

  insert_msg('Inference and Forest plots of the estimates')

  ## inference, confidence intervals are calculated based on T distribution

  incov_neuro$inference <- incov_neuro$models %>%
    map(summary,
        'inference',
        ci_method = 'distribution') %>%
    map2(., incov_neuro$variable_lexicon,
         ~mutate(.x,
                 var_lab = exchange(variable, dict = .y)))

  ## Forest plots

  incov_neuro$forest_plots <-
    list(x = incov_neuro$inference,
         plot_title = c('5-HT', 'DA sulfate'),
         x_lab = paste0(c('5-HT', 'DA sulfate'),
                        ', Z-score')) %>%
    pmap(plot_forest,
         variable = 'var_lab',
         hide_baseline = TRUE,
         estimate_size = 2.5,
         estimate_hjust = 0.5,
         cust_theme = globals$common_theme,
         plot_subtitle = paste('n =', nrow(incov_neuro$analysis_tbl)))

# Appending the Forest plots with the baselines -------

  insert_msg('Adding the baseline')

  incov_neuro$intercepts <-
    incov_neuro$inference[names(incov_neuro$forest_plots)] %>%
    map_dfr(filter,
            level == 'baseline') %>%
    mutate(baseline_lab = paste0('baseline: ', signif(estimate, 2))) %>%
    .$baseline_lab

  incov_neuro$forest_plots <-
    map2(incov_neuro$forest_plots,
         incov_neuro$intercepts,
         ~.x +
           labs(x = paste(.x$labels$x,
                          .y,
                          sep = ', ')))

# Caching the results ------

  insert_msg('Caching the results')

  save(incov_neuro, file = './cache/incov_neuro.RData')

# END ------

  insert_tail()
