# Comparison of the SIMMUN and INCOV cohort

  insert_head()

# container ------

  comp <- list()

# analysis globals ----

  insert_msg('Analysis globals')

  ## common variables

  comp$lexicon <-
    c('sex' = 'Sex',
      'age' = 'Age, years',
      'bmi_class' = 'Body mass',
      'cov' = 'SARS-CoV-2',
      'severity' = 'SARS-CoV-2 hospitalization') %>%
    compress(names_to = 'variable',
             values_to = 'label')

  ## analysis tables

  comp$analysis_tbl <-
    list(SIMMUN = stigma$data %>%
           filter(patient_id %in% stigma$complete_ids),
         INCOV = incov$clinic %>%
           filter(patient_id %in% incov$complete_ids) %>%
           mutate(severity = who_severity)) %>%
    map(filter, !duplicated(patient_id)) %>%
    map(select, all_of(comp$lexicon$variable))

  comp$analysis_tbl$INCOV <-
    comp$analysis_tbl$INCOV %>%
    mutate(severity = ifelse(cov == 'healthy',
                             'healthy',
                             ifelse(stri_detect(severity,
                                                regex = '(1|2)$'),
                                    'ambulatory', 'hospitalized')),
           severity = factor(severity, c('healthy', 'ambulatory', 'hospitalized')))

  comp$analysis_tbl <- comp$analysis_tbl %>%
    compress(names_to = 'cohort') %>%
    mutate(cohort = factor(cohort, c('SIMMUN', 'INCOV')))

  ## test type

  comp$lexicon$test_type <-
    comp$analysis_tbl[comp$lexicon$variable] %>%
    map_lgl(is.numeric) %>%
    ifelse('wilcoxon_r', 'cramer_v')

# Descriptive stats ----

  insert_msg('Descriptive stats')

  comp$stats <- comp$analysis_tbl %>%
    explore(variables = comp$lexicon$variable,
            split_factor = 'cohort',
            what = 'table',
            pub_styled = TRUE) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('variable', levels(comp$analysis_tbl$cohort)))

# Testing ------

  insert_msg('Testing')

  comp$test <- comp$analysis_tbl %>%
    compare_variables(variables = comp$lexicon$variable,
                      split_factor = 'cohort',
                      what = 'eff_size',
                      types = comp$lexicon$test_type,
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH')

# Result table ------

  insert_msg('Result table')

  comp$result_tbl <-
    left_join(comp$stats,
              comp$test[c('variable', 'test', 'significance', 'eff_size')],
              by = 'variable') %>%
    format_tbl(rm_complete = TRUE,
               dict = comp$lexicon) %>%
    full_rbind(tibble(variable = 'Patritipants, n',
                      SIMMUN = table(comp$analysis_tbl$cohort)[1],
                      INCOV = table(comp$analysis_tbl$cohort)[2]), .) %>%
    mutate(test = stri_replace(test,
                               fixed = ' test',
                               replacement = ''),
           test = stri_replace(test,
                               fixed = 'Chi-squared',
                               replacement = '\u03C7\u00B2')) %>%
    set_names(c('Variable', 'SIMMUN', 'INCOV',
                'Test', 'Significance',  'Effect size'))

# END -----

  insert_tail()
