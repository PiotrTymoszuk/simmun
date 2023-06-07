# manuscript tables

  insert_head()

# container -----

  tabs <- list()
  suppl_tabs <- list()

# Tables 1 and 2: characteristic of the STIGMA and INCOV collectives -------

  insert_msg('Table 1 and 2: characteristic of the STIGMA and INCOV collective')

  tabs[c('stigma_cohort',
         'incov_cohort')] <-
    cohort$feat_tables %>%
    map(filter,!stri_detect(Variable, regex = '(post-COVID)|(BMI)')) %>%
    map(~map_dfc(.x,
                 stri_replace,
                 regex = '^HADS.*\\nHADS\\+:\\s{1}',
                 replacement = '')) %>%
    map(set_names,
        c('Variable', 'Uninfected', 'SARS-CoV-2 infection',
          'Test', 'Significance', 'Effect size'))

  ## skipping the somatic symptoms

  tabs$stigma_cohort <- tabs$stigma_cohort %>%
    filter(Variable %in% c('Participants, n',
                           'Sex',
                           'Age, years',
                           'Body mass index',
                           'Physical disorder',
                           'Mental disorder',
                           'HADS anxiety score',
                           'HADS depression score',
                           'Clinically relevant signs of depression or anxiety, HADS ≥ 8',
                           'PSS-4 mental stress score',
                           'anti-RBD SARS-CoV-2, IgG, AU',
                           'COVID-19 severity'))

  ## the ready to use tables

  tabs[c('stigma_cohort',
         'incov_cohort')] <- tabs[c('stigma_cohort',
                                    'incov_cohort')] %>%
    list(x = .,
         label = c('table_1_stigma_cohort',
                   'table_2_incov_cohort'),
         ref_name = c('stigma_cohort',
                      'incov_cohort'),
         caption = c(paste('Characteristic of the local SIMMUN cohort.',
                           'Numeric variables are presented as medians with',
                           'interquartile ranges. Categorical variables are',
                           'presented as percentages and counts within',
                           'complete observation sets.'),
                     paste('Characteristic of the external INCOV cohort.',
                           'Numeric variables are presented as medians with',
                           'interquartile ranges. Categorical variables are',
                           'presented as percentages and counts within',
                           'complete observation sets.'))) %>%
    pmap(mdtable)

# Table S1: variables used in modeling, SIMMUN ------

  insert_msg('Table S1: study variables')

  ## study variables with their format and transformations

  suppl_tabs$study_vars <-
    list(response = stigma$responses,
         explanatory = stigma$expl_lexicon$variable) %>%
    map(~filter(globals$stigma_lexicon,
                variable %in% .x)) %>%
    compress(names_to = 'variable_type')

  suppl_tabs$study_vars$format <-
    stigma$analysis_tbl[suppl_tabs$study_vars$variable] %>%
    map_lgl(is.numeric) %>%
    ifelse('numeric', 'categorical')

  suppl_tabs$study_vars$levels <-
    stigma$analysis_tbl[suppl_tabs$study_vars$variable] %>%
    map(levels) %>%
    map_chr(function(x) if(is.null(x)) NA else paste(x, collapse = ', '))

  suppl_tabs$study_vars <- suppl_tabs$study_vars %>%
    transmute(`Variable type` = variable_type,
              `Variable label` = base_label,
              Format = format,
              Unit = unit,
              Transformation = stri_extract(axis_label, regex = '^(log|sqrt)'),
              Transformation = car::recode(Transformation,
                                           "'log' = 'logarithm'; 'sqrt' = 'square root'"),
              Transformation = ifelse(is.na(Transformation) & Format == 'numeric',
                                      'identity', Transformation),
              Categories = levels) %>%
    mdtable(label = 'table_s1_study_variables',
            ref_name = 'study_vars',
            caption = 'Study variables.')

# Table S2: comparison of the included and excluded participants -----

  insert_msg('Table S2: included and excluded participants, SIMMUN')

  suppl_tabs$excluded <- excl$feat_table %>%
    mdtable(label = 'table_s2_included_excluded',
            ref_name = 'excluded',
            caption = paste('Significant differences between participant',
                            'of the SIMMUN study included in the analysis',
                            'and SIMMUN participants excluded due to data',
                            'missingness.'))

# Table S3: study variables, INCOV -----

  insert_msg('Study variables, INCOV')

  suppl_tabs$incov_study_vars <-
    globals$incov_lexicon[c('variable', 'label', 'unit')] %>%
    rbind(tibble(variable = c('age', 'sex', 'bmi_class'),
                 label = c('age', 'sex', 'body mass index'),
                 unit = c('years', NA, NA))) %>%
    mutate(var_type = ifelse(variable %in% c('serotonin', 'dopamine 3-O-sulfate'),
                             'response', 'explanatory'),
           format = ifelse(is.na(unit), 'categorical', 'numeric')) %>%
    arrange(desc(var_type)) %>%
    select(var_type,
           label,
           format,
           unit) %>%
    set_names(c('Variable type', 'Variable label', 'Format', 'Unit')) %>%
    mdtable(label = 'table_s3_incov_study_vars',
            ref_name = 'incov_study_vars',
            caption = 'Study variables in the INCOV cohort.')

# Table S4: numbers of samples available per time point -----

  insert_msg('Table S4: samples per timepoint')

  suppl_tabs$samples <- cohort$samples %>%
    mdtable(label = 'table_s4_sample_number',
            ref_name = 'samples',
            caption = paste('Number of available samples and sampling',
                            'timepoints in the INCOV cohort.'))

# Table S5: comparison of the cohorts -------

  insert_msg('Table S5: comparison of the cohorts')

  suppl_tabs$cohort_comp <- comp$result_tbl %>%
    mdtable(label = 'table_s5_cohort_comparison',
            ref_name = 'cohort_comp',
            caption = paste('Comparison of the SIMMUN and INCOV cohorts.',
                            'Numeric variables are presented as medians with',
                            'interquartile ranges (IQR) and ranges.',
                            'Categorical variables are presented as percentages',
                            'and counts within the complete observation set.'))

# Table S6: modeling, SIMMUN ------

  insert_msg('Table S6: modeling SIMMUN')

  suppl_tabs$modeling_simmun <- mod_eli$inference %>%
    map_dfr(transmute,
            Response = exchange(response, globals$stigma_lexicon),
            `Explanatory variable` = ifelse(variable == 'Intercept',
                                            variable,
                                            exchange(variable,
                                                     dict = globals$stigma_lexicon)),
            Category = level,
            n = n,
            `Estimate, 95% CI` = paste0(signif(estimate, 2), ' [',
                                        signif(lower_ci, 2), ' - ',
                                        signif(upper_ci, 2), ']'),
            Significance = ifelse(p_value >= 0.05,
                                  paste0('ns (', signif(p_value, 2), ')'),
                                  paste('p =', signif(p_value, 2)))) %>%
    mdtable(label = 'table_s5_modeling_simmun',
            ref_name = 'modeling_simmun',
            caption = paste('Results of multi-paramater linear modeling',
                            'of serum concentrations of tryptophan,',
                            'kynunerine, tyrosine, and kynurenine - tryptophan,',
                            'and phenylalanine - tyrosine ratios in the',
                            'SIMMUN cohort.'))

# Table S7 - S8: univariate analysis ------

  insert_msg('Tables S7 - S8, univariable analysis')

  suppl_tabs[c('univariate_correlation',
               'univariate_complarison')] <-
    list(x = tst_results,
         label = c('table_s7_uni_correlation',
                   'table_s8_uni_comparison'),
         ref_name = c('univariate_correlation',
                      'univariate_complarison'),
         caption = c(paste('Correlation of age, stress scoring',
                           'blood neopterin levels and neutrophils - leukocyte',
                           'ratio with serum metabolite levels in the',
                           'SIMMUN cohort.'),
                     paste('Comparison of serum metabolite levels',
                           'in SIMMUN study participants',
                           'split by depression signs',
                           'and SARS-CoV-2 infection status.'))) %>%
    pmap(mdtable)

# Table S9: modeling, INCOV ------

  insert_msg('Table S9: robust linear modeling in the INCOV cohort')

  suppl_tabs$modeling_incov <- incov_neuro$inference %>%
    map_dfr(transmute,
            Response = exchange(response, incov_neuro$response_lexicon),
            `Explanatory variable` = ifelse(variable == 'Intercept',
                                            variable,
                                            exchange(variable,
                                                     dict = reduce(incov_neuro$variable_lexicons,
                                                                   rbind))),
            Category = level,
            n = n,
            `Estimate, 95% CI` = paste0(signif(estimate, 2), ' [',
                                        signif(lower_ci, 2), ' - ',
                                        signif(upper_ci, 2), ']'),
            Significance = ifelse(p_value >= 0.05,
                                  paste0('ns (', signif(p_value, 2), ')'),
                                  paste('p =', signif(p_value, 2)))) %>%
    mdtable(label = 'table_s9_modeling_simmun',
            ref_name = 'modeling_incov',
            caption = paste('Results of multi-paramater robust linear modeling',
                            'of serum concentrations of serotonin,',
                            'and dopamine sulfate in the',
                            'INCOV cohort.'))

# Table S10 - S11: time-course modeling -----

  insert_msg('Table S10 - S11: time course modeling')

  suppl_tabs[c('time_uninfected',
               'time_acute')] <-
    list(x = time_rlm$result_tbl,
         label = c('table_s11_time_modeling_incov_uninfected',
                   'table_s12_time_modeling_incov_acute'),
         ref_name = c('time_uninfected',
                      'time_acute'),
         caption = c(paste('Results of robust linear modeling',
                           'of serum levels of inflammatory',
                           'cytokines, tryptophan, tyrosine and',
                           'their metabolites as a function of',
                           'SARS-CoV-2 infection status in the INCOV cohort.',
                           'The uninfected subset served as a baseline.'),
                     paste('Results of robust linear modeling',
                           'of serum levels of inflammatory',
                           'cytokines, tryptophan, tyrosine and',
                           'their metabolites as a function of',
                           'SARS-CoV-2 infection status in the INCOV cohort.',
                           'The acute SARS-CoV-2 infection subset',
                           'served as a baseline.'))) %>%
    pmap(mdtable)

# END ------

  insert_tail()
