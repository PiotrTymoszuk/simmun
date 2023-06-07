# Age,metabolites and cytokines for particular infection timepoints

  insert_head()

# container ------

  age_incov <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  # analysis globals ------

  insert_msg('Analysis globals')

  ## variable lexicon

  age_incov$lexicon <- globals$incov_lexicon %>%
    mutate(label_long = variable,
           variable = make.names(variable))

  ## analysis tables

  age_incov$analysis_tbl <- incov$analysis_tbl %>%
    select(patient_id,
           timepoint,
           age,
           any_of(set_names(age_incov$lexicon$label_long,
                            age_incov$lexicon$variable))) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"),
           timepoint = factor(timepoint,
                              c('uninfected', 'acute', 'sub-acute', 'recovery'))) %>%
    blast(timepoint) %>%
    map(~filter(.x, complete.cases(.x)))

# Correralation analysis -----

  insert_msg('Correlation analysis')

  age_incov$test <- age_incov$analysis_tbl %>%
    map(function(data) map(age_incov$lexicon$variable, ~c('age', .x)) %>%
          map_dfr(~correlate_variables(data,
                                       variables = .x,
                                       what = 'correlation',
                                       type = 'spearman',
                                       ci = TRUE,
                                       pub_styled = FALSE))) %>%
    compress(names_to = 'timepoint') %>%
    re_adjust

# Bubble plot -----

  insert_msg('Bubble plot')

  age_incov$plot <- age_incov$test %>%
    mutate(variable1 = timepoint) %>%
    corr_buble(plot_title = 'Metabolites, cytokines and age',
               signif_only = FALSE) +
    labs(plot_subtitle = 'Spearman correlation with age') +
    scale_y_continuous(labels = function(x) exchange(x, dict = age_incov$lexicon))

# END ------

  insert_tail()
