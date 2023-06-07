# Testing for significant association of single factors for the TRP/KYN system

  insert_head()

# container -----

  tst_trp <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variables

  tst_trp$lexicon <- globals$stigma_lexicon %>%
    filter(variable %in% test_globals$variables$trp)

  tst_trp$responses <- stigma$responses[c('trp', 'kyn', 'kyn_trp')]

  ## analysis tables with non-normalized variable values

  tst_trp$analysis_tbl <- stigma$data %>%
    mutate(infection = car::recode(cov,
                                   "'healthy' = 'uninfected';
                                    'SARS-CoV-2' = 'SARS-CoV-2'"),
           infection = factor(infection,
                              c('uninfected', 'SARS-CoV-2'))) %>%
    select(patient_id,
           all_of(unname(tst_trp$responses)),
           all_of(tst_trp$lexicon$variable)) %>%
    filter(patient_id %in% stigma$complete_ids)

  ## analysis types

  tst_trp$lexicon$analysis_type <-
    tst_trp$analysis_tbl[c(tst_trp$lexicon$variable)] %>%
    map_lgl(is.numeric) %>%
    ifelse('correlation', 'comparison')

  ## variable types

  tst_trp$var_types <- tst_trp$lexicon %>%
    blast(analysis_type) %>%
    map(~set_names(.x$variable,
                   .x$variable))

# Comparison ------

  insert_msg('Comparison')

  ## descriptive stats

  tst_trp$comparison$stats <- tst_trp$var_types$comparison %>%
    map(~explore(tst_trp$analysis_tbl,
                 split_factor = .x,
                 variables = tst_trp$responses,
                 what = 'table',
                 pub_styled = TRUE)) %>%
    map(compress, names_to = 'level') %>%
    compress(names_to = 'split_factor')

  ## testing for differences: T test

  tst_trp$comparison$test <- tst_trp$var_types$comparison %>%
    map(~compare_variables(tst_trp$analysis_tbl,
                           variables = tst_trp$responses,
                           split_factor = .x,
                           what = 'eff_size',
                           types = 'cohen_d',
                           ci = FALSE,
                           pub_styled = TRUE)) %>%
    compress(names_to = 'split_factor') %>%
    mutate(comp_id = paste(split_factor, variable, sep = '|'))

# Correlation ------

  insert_msg('Correlation')

  ## variable pairs

  tst_trp$correlation$pairs <- tst_trp$var_types$correlation %>%
    map(function(x) tst_trp$responses %>%
          map(~c(x, .x))) %>%
    unlist(recursive = FALSE)

  tst_trp$correlation$pairs <- tst_trp$correlation$pairs %>%
    map_chr(~paste(.x[1], .x[2], sep = '|')) %>%
    set_names(tst_trp$correlation$pairs, .)

  ## Pearson's correlations

  tst_trp$correlation$test <- tst_trp$correlation$pairs %>%
    map_dfr(~correlate_variables(tst_trp$analysis_tbl,
                                 variables = .x,
                                 what = 'correlation',
                                 type = 'pearson',
                                 ci = TRUE,
                                 pub_styled = FALSE)) %>%
    re_adjust %>%
    mutate(comp_id = paste(variable1, variable2, sep = '|'), )

# Multiple testing adjustments ------

  insert_msg('FDR adjustment')

  ## global adjustment for multiple comparison with the FDR method

  tst_trp$fdr <- list(comparison = tst_trp$comparison$test,
                      correlation = tst_trp$correlation$test) %>%
    map(~.x[c('comp_id', 'p_value')]) %>%
    compress(names_to = 'analysis_type') %>%
    re_adjust %>%
    mutate(analysis_type = factor(analysis_type,
                                  c('comparison', 'correlation')))

  ## merging with the genuine testing results

  tst_trp$comparison$test <-
    left_join(tst_trp$comparison$test %>%
                select(-p_value, -p_adjusted, -significance),
              tst_trp$fdr,
              by = 'comp_id') %>%
    mutate(eff_size = stri_replace(eff_size, fixed = '-', replacement = ''),
           plot_cap = paste(eff_size, significance, sep = ', '))

  tst_trp$correlation$test <-
    left_join(tst_trp$correlation$test %>%
                select(-p_value, -p_adjusted, -significance),
              tst_trp$fdr,
              by = 'comp_id') %>%
    mutate(est_lab = paste0(signif(estimate, 2), ' [',
                            signif(lower_ci, 2), ' - ',
                            signif(upper_ci, 2), ']'),
           plot_cap = paste(est_lab, significance, sep = ', '))

# Visualization of the comparisons -----

  insert_msg('Visualiztion: comparisons')

  tst_trp$comparison$plots <-
    list(variable = tst_trp$comparison$test$variable,
         split_factor = tst_trp$comparison$test$split_factor,
         plot_title = exchange(tst_trp$comparison$test$variable,
                               dict = globals$stigma_lexicon,
                               value = 'base_label'),
         plot_subtitle = tst_trp$comparison$test$plot_cap,
         y_lab = exchange(tst_trp$comparison$test$variable,
                          dict = globals$stigma_lexicon,
                          value = 'axis_label'),
         x_lab = exchange(tst_trp$comparison$test$split_factor,
                          dict = tst_trp$lexicon)) %>%
    pmap(plot_variable,
         tst_trp$analysis_tbl,
         type = 'box',
         cust_theme = globals$common_theme,
         x_n_labs = TRUE) %>%
    map(~.x +
          scale_fill_brewer(palette = 'Reds') +
          theme(legend.position = 'none')) %>%
    set_names(tst_trp$comparison$test$comp_id)

# visualization of the correlations -----

  insert_msg('Visualization of the correlations')

  tst_trp$correlation$plots <-
    list(variables = map2(tst_trp$correlation$test$variable1,
                          tst_trp$correlation$test$variable2,
                          c),
         plot_title = map2(exchange(tst_trp$correlation$test$variable1,
                                    dict = globals$stigma_lexicon,
                                    value = 'base_label'),
                           exchange(tst_trp$correlation$test$variable2,
                                    dict = globals$stigma_lexicon,
                                    value = 'base_label'),
                           paste, sep = ' and '),
         plot_subtitle = tst_trp$correlation$test$plot_cap,
         x_lab = exchange(tst_trp$correlation$test$variable1,
                          dict = globals$stigma_lexicon,
                          value = 'axis_label'),
         y_lab = exchange(tst_trp$correlation$test$variable2,
                          dict = globals$stigma_lexicon,
                          value = 'axis_label')) %>%
    pmap(plot_correlation,
         tst_trp$analysis_tbl,
         type = 'correlation',
         point_color = 'indianred3',
         cust_theme = globals$common_theme) %>%
    map(~.x +
          theme(plot.tag = element_blank())) %>%
    set_names(tst_trp$correlation$test$comp_id)

# Correlograms ------

  insert_msg('Correlaograms')

  tst_trp$correlation$correlogram <- tst_trp$correlation$test %>%
    corr_buble(signif_only = FALSE,
               plot_title = 'TRP - KYN pathway',
               rotate_x_labs = FALSE,
               fill_title = 'r') +
    scale_x_discrete(labels = function(x) exchange(x, dict = globals$stigma_lexicon) %>%
                       stri_replace(fixed = ' score', replacement = '')) +
    scale_y_discrete(labels = function(x) exchange(x, dict = globals$stigma_lexicon) %>%
                       stri_replace(fixed = ' score', replacement = ''),
                     limits = tst_trp$responses)

# END ------

  insert_tail()
