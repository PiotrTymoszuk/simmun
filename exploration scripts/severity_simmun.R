# Differences in severity of COVID-19 attributed to age and gender

  insert_head()

# container ----

  sev_simmun <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  sev_simmun$sev_colors <- c(healthy = 'steelblue',
                             ambulatory = 'coral3',
                             hospitalized = 'coral4')

# age ------

  insert_msg('Age nd severity')

  ## stats

  sev_simmun$age$stats <- stigma$data %>%
    explore(split_factor = 'severity',
            variables = 'age',
            what = 'table',
            pub_styled = TRUE) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('variable', levels(stigma$data$severity)))

  ## testing for differences in age, Kruskal-Wallis test

  sev_simmun$age$test <- stigma$data %>%
    compare_variables(variables = 'age',
                      split_factor = 'severity',
                      what = 'eff_size',
                      types = 'kruskal_etasq',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

  ## plot

  sev_simmun$age$plot <- stigma$data %>%
    plot_variable(variable = 'age',
                  split_factor = 'severity',
                  type = 'box',
                  point_hjitter = 0,
                  cust_theme = globals$common_theme,
                  plot_title = 'Age and SARS-CoV-2 infection severity, SIMMUN',
                  plot_subtitle = sev_simmun$age$test$plot_cap,
                  x_n_labs = TRUE) +
    scale_fill_manual(values = sev_simmun$sev_colors)

# gender -----

  insert_msg('Severity and gender')

  ## stats

  sev_simmun$sex$stats <-  stigma$data %>%
    explore(split_factor = 'sex',
            variables = 'severity',
            what = 'table',
            pub_styled = TRUE) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('variable', levels(stigma$data$sex)))

  ## testing: chi-squared test with Cramer's V effect size

  sev_simmun$sex$test <- stigma$data %>%
    compare_variables(variables = 'severity',
                      split_factor = 'sex',
                      what = 'eff_size',
                      types = 'cramer_v',
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE) %>%
    mutate(plot_cap = paste(eff_size, significance, sep = ', '))

  ## plot

  sev_simmun$sex$plot <- stigma$data %>%
    plot_variable(variable = 'severity',
                  split_factor = 'sex',
                  type = 'stack',
                  scale = 'percent',
                  cust_theme = globals$common_theme,
                  plot_title = 'Gender and SARS-CoV-2 infection severity',
                  plot_subtitle = sev_simmun$sex$test$plot_cap,
                  x_n_labs = TRUE) +
    scale_fill_manual(values = sev_simmun$sev_colors)

# END ------

  insert_tail()
