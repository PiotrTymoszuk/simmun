# Normality of distribution and equality of variances
#
# Split by the the study group (SARS-CoV2 status/HADS, STIGMA)
# and timepoint (INCOV)

  insert_head()

# container -----

  distr <- list()

# globals: analysis tables -------

  insert_msg('Analysis tables')

  ## variables, in the SIMMUN/STIGMA cohort, CRP and IL6
  ## are excluded from testing, they are handled as factors in modeling

  distr$variables <- map(globals[c('stigma_lexicon',
                                   'incov_lexicon')],
                         ~.x$variable) %>%
    set_names(c('stigma', 'incov'))

  distr$variables$stigma <-
    distr$variables$stigma[!distr$variables$stigma %in% c('crp', 'il6', 'infection')]

  distr$variables$stigma <-
    distr$variables$stigma[map_lgl(stigma$data[distr$variables$stigma],
                                   is.numeric)]

  ## analysis tables for both cohorts

  distr$analysis_tbl$stigma <- stigma$data

  distr$analysis_tbl$incov <- incov[c('proteome',
                                      'metabolome')] %>%
    map(select,
        sample_id,
        patient_id,
        timepoint,
        any_of(distr$variables$incov))

  ## reducing the INCOV dataset to a single data frame

  distr$analysis_tbl$incov <- distr$analysis_tbl$incov %>%
    reduce(left_join,
           by = c('sample_id', 'patient_id', 'timepoint'))

# Normality: Shapiro-Wilk tests ------

  insert_msg('Normality: Shapiro-Wilk test')

  distr$shapiro_test <-
    list(data = distr$analysis_tbl,
         variables = distr$variables,
         split_factor = c('study_group', 'timepoint')) %>%
    pmap(explore,
         what = 'normality',
         pub_styled = FALSE)

  ## indicating the best transformation of the STIGMA variables
  ## in terms of normality

  distr$stigma_normality_trans <-
    distr$shapiro_test$stigma %>%
    map(mutate,
        variable = factor(variable, distr$variables$stigma),
        base_variable = stri_replace(variable,
                                     regex = '(log_)|(sqrt_)',
                                     replacement = ''),
        transformation = ifelse(stri_detect(variable, fixed = 'log'),
                                'log',
                                ifelse(stri_detect(variable, fixed = 'sqrt'),
                                       'sqrt', 'identity'))) %>%
    map(ddply, 'base_variable', filter, stat == max(stat)) %>%
    map(select, variable, base_variable, transformation, stat) %>%
    map(as_tibble)

  distr$stigma_normality_best <- distr$stigma_normality_trans %>%
    compress(names_to = 'study_group') %>%
    ddply('base_variable', filter, stat == max(stat)) %>%
    as_tibble

# EOV: Levene test --------

  insert_msg('Homogeneity of variances')

  ## STIGMA: no hard violations detected

  distr$levene_test <-
    list(data = distr$analysis_tbl,
         var = distr$variables,
         splt = c('study_group', 'timepoint')) %>%
    pmap(function(data, var, splt) compare_variables(data,
                                                     variables = var,
                                                     split_factor = splt,
                                                     what = 'variance',
                                                     pub_styled = TRUE))

# END -------

  insert_tail()
