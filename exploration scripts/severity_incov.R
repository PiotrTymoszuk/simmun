# Effects of age an gender on bi_severity of COVID-19

  insert_head()

# container -----

  sev_incov <- list()

# analysis globals -----

  insert_msg('Analysis globals')

  sev_incov$analysis_tbl <- incov$clinic %>%
    filter(patient_id %in% unique(incov$complete_ids)) %>%
    filter(!duplicated(patient_id), !is.na(bi_severity)) %>%
    mutate(bi_severity = droplevels(bi_severity)) %>%
    select(patient_id, age, sex, bi_severity)

  sev_incov$sev_colors <-
    c(healthy = 'steelblue',
      mild = 'darkolivegreen4',
      moderate = 'coral3',
      severe = 'coral4',
      critical = 'gray40',
      `mild/moderate` = 'steelblue',
      `severe/critical` = 'coral3')

# Differences in age ------

  insert_msg('Age')

  ## stats

  sev_incov$age$stats <- sev_incov$analysis_tbl %>%
    explore(split_factor = 'bi_severity',
            variables = 'age',
            what = 'table',
            pub_styled = TRUE)

  ## Kruskal-Wallis test

  sev_incov$age$test <- sev_incov$analysis_tbl %>%
    compare_variables(variables = 'age',
                      split_factor = 'bi_severity',
                      what = 'eff_size',
                      types = 'kruskal_etasq',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

  ## plot

  sev_incov$age$plot <- sev_incov$analysis_tbl %>%
    plot_variable(variable = 'age',
                  split_factor = 'bi_severity',
                  type = 'box',
                  point_hjitter = 0,
                  cust_theme = globals$common_theme,
                  plot_title = 'Age and SARS-CoV-2 infection bi_severity',
                  plot_subtitle = sev_incov$age$test$plot_cap,
                  x_n_labs = TRUE) +
    scale_fill_manual(values = sev_incov$sev_colors)

# Gender and bi_severity ------

  insert_msg('Gender and bi_severity')

  ## descriptive stats

  sev_incov$sex$stats <- sev_incov$analysis_tbl %>%
    explore(split_factor = 'sex',
            variables = 'bi_severity',
            what = 'table',
            pub_styled = TRUE) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('variable', levels(sev_incov$analysis_tbl$sex)))

  ## Chi-squared test

  sev_incov$sex$test <- sev_incov$analysis_tbl %>%
    compare_variables(variables = 'bi_severity',
                      split_factor = 'sex',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

  ## plot

  sev_incov$sex$plot <- sev_incov$analysis_tbl %>%
    plot_variable(variable = 'bi_severity',
                  split_factor = 'sex',
                  scale = 'percent',
                  type = 'stack',
                  cust_theme = globals$common_theme,
                  plot_title = 'Gender and SASR-CoV-2 bi_severity, INCOV',
                  plot_subtitle = sev_incov$sex$test$plot_cap,
                  x_n_labs = TRUE) +
    scale_fill_manual(values = sev_incov$sev_colors)

# END ------

  insert_tail()
