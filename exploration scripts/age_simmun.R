# Effects of age on metabolites and markers of inflammation
# correlation analysis (Pearson)

  insert_head()

# container ------

  age_simmun <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## variables

  age_simmun$variables <- c(neo = 'log_neo',
                               nlr = 'log_nlr',
                               stigma$responses)

  ## analysis table

  age_simmun$analysis_tbl <-
    stigma$data[c('patient_id', 'age', age_simmun$variables)] %>%
    filter(complete.cases(.))

# Correlation analysis ------

  insert_msg('Correlation analysis')

  age_simmun$test <- age_simmun$variables %>%
    map_dfr(~correlate_variables(age_simmun$analysis_tbl,
                                 variables = c('age', .x),
                                 what = 'correlation',
                                 type = 'pearson',
                                 ci = TRUE,
                                 pub_styled = TRUE)) %>%
    re_adjust %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '),
           plot_cap = paste(plot_cap, n, sep = ', n = '))

# Plotting -------

  insert_msg('Plotting, scatter plots')

  age_simmun$plots <-
    list(variables = map(age_simmun$variables, ~c('age', .x)),
         plot_title = exchange(names(age_simmun$variables),
                               globals$stigma_lexicon),
         plot_subtitle = age_simmun$test$plot_cap,
         y_lab = exchange(age_simmun$variables,
                          globals$stigma_lexicon,
                          value = 'axis_label')) %>%
    pmap(plot_correlation,
         age_simmun$analysis_tbl,
         type = 'correlation',
         cust_theme = globals$common_theme,
         x_lab = 'age, years')

# END ------

  insert_tail()
