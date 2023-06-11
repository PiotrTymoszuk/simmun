# Comparison of the included and excluded SIMMUN/STIGMA participants

  insert_head()

# container ----

  excl <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variable lexicon

  excl$var_lexicon <-
    c(sex = 'Sex',
      age = 'Age, years',
      bmi_class = 'Body mass class',
      somatic_comorb = 'Physical illness',
      psych_comorb = 'Mental illness',
      hads_anx_score = 'HADS anxiety score',
      hads_dpr_score = 'HADS depression score',
      hads_signs = 'Depression or anxiety signs, HADS \u2265 8',
      pss_stress_score = 'PSS-4 stress score',
      cov = 'SARS-CoV-2 infection',
      anti_rbd = 'anti-RBD SARS-CoV-2, IgG, AU',
      severity = 'COVID-19 severity',
      trp = 'TRP, µmol/L',
      kyn = 'KYN, µmol/L',
      kyn_trp = 'KYN/TRP',
      tyr = 'TYR, µmol/L',
      phe = 'PHE, µmol/L',
      phe_tyr = 'PHE/TYR',
      neo = 'NEO, nmol/L',
      nlr = 'NLR') %>%
    compress(names_to = 'variable',
             values_to = 'label')

  excl$var_lexicon$test_type <- stigma$data[excl$var_lexicon$variable] %>%
    map_lgl(is.numeric) %>%
    ifelse('wilcoxon_r', 'cramer_v')

  ## analysis table

  excl$analysis_tbl <- stigma$data %>%
    mutate(analysis_status = ifelse(patient_id %in% stigma$complete_ids,
                                    'analyzed', 'excluded'),
           analysis_status = factor(analysis_status,
                                    c('analyzed', 'excluded'))) %>%
    select(patient_id,
           analysis_status,
           all_of(excl$var_lexicon$variable)) %>%
    mutate(cov = car::recode(cov,
                             "'healthy' = 'uninfected'"),
           cov = factor(cov, c('uninfected', 'SARS-CoV-2')),
           severity = car::recode(severity,
                                  "'healthy' = 'uninfected'"),
           severity = factor(severity,
                             c('uninfected', 'ambulatory', 'hospitalized')))

# Descriptive stats -----

  insert_msg('Descriptive stats')

  excl$stats <- excl$analysis_tbl %>%
    explore(split_factor = 'analysis_status',
            variables = excl$var_lexicon$variable,
            what = 'table',
            pub_styled = TRUE) %>%
    reduce(left_join, by = 'variable') %>%
    set_names(c('variable', levels(excl$analysis_tbl$analysis_status)))

# Testing -------

  insert_msg('Testing')

  excl$test <- excl$analysis_tbl %>%
    compare_variables(variables = excl$var_lexicon$variable,
                      split_factor = 'analysis_status',
                      what = 'eff_size',
                      types = excl$var_lexicon$test_type,
                      exact = FALSE,
                      ci = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH')

# A table with the significant results ------

  insert_msg('Result table')

  excl$feat_table <-
    left_join(excl$stats,
              excl$test[c('variable', 'test', 'significance', 'eff_size')],
              by = 'variable') %>%
    filter(!stri_detect(significance, regex = '^ns')) %>%
    format_tbl(rm_complete = FALSE,
               dict = excl$var_lexicon) %>%
    mutate(test = stri_replace(test,
                               fixed = ' test',
                               replacement = ''),
           test = stri_replace(test,
                               fixed = 'Chi-squared',
                               replacement = '\u03C7\u00B2')) %>%
    set_names(c('Variable', 'Analyzed', 'Excluded',
                'Test', 'Significance', 'Effect size'))

# END ------

  insert_tail()
