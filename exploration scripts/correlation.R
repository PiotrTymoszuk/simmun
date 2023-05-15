# Inter-correlations of each of the outcomes, inflammatory markers
# and psychometric scores

  insert_head()

# container ------

  corr <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  corr$variables <-
    list(outcomes = c('trp', 'log_kyn', 'log_kyn_trp',
                      'log_phe', 'log_tyr', 'sqrt_phe_tyr'),
         inflammation = c('log_neo', 'log_nlr'),
         psychometry = c('pss_stress_score',
                         'hads_anx_score',
                         'hads_dpr_score'))

  corr$pairs <- corr$variables %>%
    map(~combn(.x, m = 2, simplify = FALSE))

# Correlation -----

  insert_msg('Correlation')

  corr$test <-
    list(x = corr$pairs,
         y = c('pearson', 'pearson', 'spearman')) %>%
    pmap(function(x, y) x %>%
           map_dfr(~correlate_variables(stigma$analysis_tbl,
                                        variables = .x,
                                        what = 'correlation',
                                        type = y,
                                        pub_styled = FALSE))) %>%
    compress(names_to = 'task') %>%
    mutate(task = factor(task, names(corr$variables))) %>%
    re_adjust %>%
    mutate(significant = ifelse(p_adjusted < 0.05,
                                'yes', 'no'),
           fontface = ifelse(significant == 'yes',
                             'bold', 'plain')) %>%
    blast(task)

# Correlograms ------

  insert_msg('Correlograms')

  corr$plots <-
    list(x = corr$test,
         y = c('Metabolites',
               'Inflammatory markers',
               'Psychometry'),
         z = corr$variables,
         v = c('Pearson correlation',
               'Pearson correlation',
               'Spearman correlation')) %>%
    pmap(function(x, y, z, v) x %>%
           ggplot(aes(x = variable1,
                      y = variable2,
                      fill = estimate,
                      size = abs(estimate)))  +
           geom_point(shape = 21) +
           geom_text(aes(label = signif(estimate, 2),
                         fontface = fontface,
                         alpha = significant),
                     size = 2.75,
                     vjust = -1.6,
                     hjust = 0.5) +
           scale_alpha_manual(values = c(yes = 1,
                                         no = 0.25)) +
           scale_fill_gradient2(low = 'steelblue',
                                mid = 'white',
                                high = 'firebrick',
                                midpoint = 0,
                                limits = c(-1, 1),
                                name = 'Corr. coefficient') +
           scale_size_area(max_size = 4.5,
                           limits = c(0, 1),
                           name = 'abs(corr. coefficient)') +
           scale_x_discrete(limits = z,
                            labels = function(x) exchange(x, dict = globals$stigma_lexicon)) +
           scale_y_discrete(limits = z,
                            labels = function(x) exchange(x, dict = globals$stigma_lexicon)) +
           guides(alpha = 'none') +
           globals$common_theme +
           theme(axis.title = element_blank()) +
           labs(title = y,
                subtitle = v))

# END ------

  insert_tail()
