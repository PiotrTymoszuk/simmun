# Characteristic of the local STIGMA and INCOV study cohorts

  insert_head()

# container ------

  cohort <- list()

# Globals -------

  insert_msg('Globals')

  ## variable lexicon

  cohort$var_lexicon <-
    c(sex = 'Sex',
      age = 'Age, years',
      bmi = 'BMI, kg/m\u00B2',
      bmi_class = 'Body mass class',
      ethnics = 'Ethnics',
      somatic_comorb = 'Physical disorder',
      psych_comorb = 'Psychiatric disorder',
      hads_anx_score = 'HADS anxiety score',
      hads_dpr_score = 'HADS depression score',
      hads_signs = 'Depression or anxiety signs, HADS \u2265 8',
      pss_stress_score = 'PSS-4 stress score',
      cov = 'Infection',
      anti_rbd = 'anti-RBD SARS-CoV-2, IgG, AU',
      severity = 'COVID-19 severity') %>%
    compress(names_to = 'variable',
             values_to = 'label')

  ## analysis tables: STIGMA

  cohort$analysis_tbl$stigma <- stigma$data %>%
    select(cov,
           patient_id,
           any_of(cohort$var_lexicon$variable)) %>%
    filter(patient_id %in% stigma$complete_ids)

  ## analysis tables: INCOV

  cohort$analysis_tbl$incov <- incov$clinic %>%
    filter(patient_id %in% incov$complete_ids) %>%
    filter(!duplicated(patient_id)) %>%
    mutate(long_cov = ifelse(long_cov == 'healthy',
                             'no', as.character(long_cov)),
           long_cov = factor(long_cov),
           psy_long_cov = ifelse(psy_long_cov %in% c('healthy',
                                                     'recovered'),
                                 'no', as.character(psy_long_cov)),
           psy_long_cov = factor(psy_long_cov),
           neuro_long_cov = ifelse(neuro_long_cov %in% c('healthy',
                                                         'recovered'),
                                   'no', as.character(neuro_long_cov)),
           neuro_long_cov = factor(neuro_long_cov),
           gastro_long_cov = ifelse(gastro_long_cov %in% c('healthy',
                                                           'recovered'),
                                    'no', as.character(gastro_long_cov)),
           gastro_long_cov = factor(gastro_long_cov),
           cogn_long_cov = ifelse(cogn_long_cov %in% c('healthy',
                                                       'recovered'),
                                  'no', as.character(cogn_long_cov)),
           cogn_long_cov = factor(cogn_long_cov),
           fatigue_long_cov = ifelse(fatigue_long_cov %in% c('healthy',
                                                             'recovered'),
                                     'no', as.character(fatigue_long_cov)),
           fatigue_long_cov = factor(fatigue_long_cov),
           resp_long_cov = ifelse(resp_long_cov %in% c('healthy',
                                                       'recovered'),
                                  'no', as.character(resp_long_cov)),
           resp_long_cov = factor(resp_long_cov))

  ## cohort-specific variable sets

  cohort$variables <- cohort$analysis_tbl %>%
    map(names) %>%
    map(~cohort$var_lexicon$variable[cohort$var_lexicon$variable %in% .x])

  ## test type

  cohort$test_type <-
    map2(cohort$variables,
         cohort$analysis_tbl,
         function(var, data) data[var] %>%
           map_chr(~ifelse(is.numeric(.x),
                           'wilcoxon_r',
                           'cramer_v')))

# Clinical stats -------

  insert_msg('Clinical characteristic')

  cohort$clin_stats <-
    map2(cohort$analysis_tbl,
         cohort$variables,
         ~explore(.x,
                  variables = .y,
                  split_factor = 'cov',
                  what = 'table',
                  pub_styled = TRUE)) %>%
    map(reduce, left_join, by = 'variable') %>%
    map(set_names, c('variable', 'Healthy', 'SARS-CoV-2'))

# Testing for differences between the healthy and SARS-CoV-2 ----

  insert_msg('Testing for differences between healthy and infected')

  cohort$clin_test <-
    list(data = cohort$analysis_tbl,
         var = cohort$variables,
         test = cohort$test_type) %>%
    pmap(function(data, var, test) compare_variables(data,
                                                     variables = var,
                                                     split_factor = 'cov',
                                                     what = 'eff_size',
                                                     types = test,
                                                     exact = FALSE,
                                                     ci = FALSE,
                                                     pub_styled = TRUE,
                                                     adj_method = 'BH'))

# Common tables with with clinical features and testing results ------

  insert_msg('A common table with clinical features')

  cohort$feat_tables <-
    map2(cohort$clin_stats,
         map(cohort$clin_test,
             ~.x[c('variable', 'test', 'significance', 'eff_size')]),
         left_join, by = 'variable') %>%
    map(mutate,
        Healthy = ifelse(stri_detect(variable,
                                     regex = '(long_cov)|(severity)'),
                         NA, Healthy),
        test = ifelse(stri_detect(variable,
                                  regex = '(long_cov)|(severity)'),
                      NA, test),
        significance = ifelse(stri_detect(variable,
                                          regex = '(long_cov)|(severity)'),
                              NA, significance),
        `SARS-CoV-2` = stri_replace(`SARS-CoV-2`,
                                    regex = '^healthy.*\\n',
                                    replacement = ''),
        test = stri_replace(test,
                            fixed = 'test',
                            replacement = ''),
        test = stri_replace(test,
                            fixed = 'Chi-squared',
                            replacement = '\u03C7\u00B2')) %>%
    map(~map_dfc(.x,
                 stri_replace_all,
                 regex = 'PSS\\-.*\\nPSS\\+:\\s{1}',
                 replacement = ''))

  cohort$feat_tables <- cohort$feat_tables %>%
    map(format_tbl,
         dict = cohort$var_lexicon,
        rm_complete = TRUE) %>%
    map(mutate,
        eff_size = ifelse(is.na(significance), NA, eff_size)) %>%
    map(set_names,
        c('Variable', 'Healthy', 'SARS-CoV-2',
          'Test', 'Significance', 'Effect size'))

  for(i in names(cohort$feat_tables)) {

    cohort$feat_tables[[i]] <- cohort$feat_tables[[i]] %>%
      full_rbind(tibble(Variable = 'Participants, n',
                        Healthy = table(cohort$analysis_tbl[[i]]$cov)[1],
                        `SARS-CoV-2` = table(cohort$analysis_tbl[[i]]$cov)[2]),
                 .)

  }

# Sample numbers -------

  insert_msg('Numbers of samples')

  ## sampling time points

  cohort$sampling_time <- incov$analysis_tbl %>%
    group_by(timepoint) %>%
    summarize(median = median(time_po),
              lower_q = quantile(time_po, 0.25, na.rm = TRUE),
              upper_q = quantile(time_po, 0.75, na.rm = TRUE)) %>%
    mutate(days_pi = paste0(round(median),
                            ' [', round(lower_q),
                            ' - ', round(upper_q),
                            ']'),
           days_pi = ifelse(timepoint == 'healthy', NA, days_pi))

  ## sample numbers

  cohort$samples <- incov$analysis_tbl %>%
    count(timepoint)

  ## common table

  cohort$samples <- left_join(cohort$sampling_time[c('timepoint', 'days_pi')],
                              cohort$samples,
                              by = 'timepoint') %>%
    set_names(c('Time point',
                'Days post infection',
                'Sample number'))
# END -----

  rm(i)

  insert_tail()
