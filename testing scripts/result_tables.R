# Tables with results of the univariate analysis

  insert_head()

# container ------

  tst_results <- list()

# correlation results -----

  insert_msg('Correlation result table')

  tst_results$correlation <- list(trp = tst_trp,
                                  tyr = tst_tyr) %>%
    map(~.x$correlation$test) %>%
    compress(names_to = 'metabolic_subsystem') %>%
    transmute(Metabolite = exchange(variable2,
                                    dict = globals$stigma_lexicon),
              `Explanatory variable` = exchange(variable1,
                                                dict = globals$stigma_lexicon),
              n = n,
              `Correlation coefficient, 95% CI` = est_lab,
              Significance = significance)

# Comparison results -----

  insert_msg('Comparison results')

  tst_results$comparison <- list(trp = tst_trp,
                                 tyr = tst_tyr) %>%
    map(~left_join(.x$comparison$stats[c('variable', 'split_factor',
                                         'level', 'statistic')],
                   .x$comparison$test[c('variable', 'split_factor',
                                        'significance', 'eff_size')],
                   by = c('variable', 'split_factor'))) %>%
    compress(names_to = 'metabolic_subsystem') %>%
    format_tbl(dict = globals$stigma_lexicon) %>%
    mutate(split_factor = exchange(split_factor, dict = globals$stigma_lexicon),
           n = stri_extract(statistic, regex = '\\d+$'),
           statistic = stri_replace(statistic, regex = '\\ncomplete.*$', replacement = ''),
           variable = factor(variable,
                             exchange(c(tst_trp$responses, tst_tyr$responses),
                                      dict = globals$stigma_lexicon))) %>%
    arrange(variable, split_factor) %>%
    select(variable,
           split_factor,
           level,
           n,
           statistic,
           significance,
           eff_size) %>%
    set_names(c('Metabolite',
                'Explanatory variable',
                'Category',
                'N',
                'Median, IQR, range',
                'Significance',
                'Effect size'))

# END ------

  insert_tail()
