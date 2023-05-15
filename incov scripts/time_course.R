# Comparison of time_course and the timepoints
# Comparison via Kruskal-Wallis test, since
# multiple variables are non-normally distributed
# and violations of EOV were detected by the Levene test

  insert_head()

# container -----

  time_course <- list()

# analysis table and variables -----

  insert_msg('Globals')

  time_course$analysis_tbl <- incov$analysis_tbl %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"),
           timepoint = factor(timepoint,
                              c('uninfected', 'acute', 'sub-acute', 'recovery')))

# Descriptive statistics --------

  insert_msg('Descriptive stats')

  time_course$stats <- time_course$analysis_tbl %>%
    explore(variables = globals$incov_lexicon$variable,
            split_factor = 'timepoint',
            what = 'table',
            pub_styled = TRUE)

# Post-hoc Mann-Whitney tests ------

  insert_msg('Post-hoc Mann-Whitney tests')

  ## the variable name 'dopamine 3-O-sulfate' is not compatible with the
  ## compare_means() function and need to be substituted by 'DA'

  time_course$post_hoc_formulas <-
    c(globals$incov_lexicon$variable[globals$incov_lexicon$variable != 'dopamine 3-O-sulfate'],
      'DA') %>%
    paste('~ timepoint') %>%
    map(as.formula) %>%
    set_names(globals$incov_lexicon$variable)

  time_course$post_hoc_test <- time_course$post_hoc_formulas %>%
    map(compare_means,
        data = time_course$analysis_tbl %>%
          mutate(DA = `dopamine 3-O-sulfate`),
        method = 'wilcox.test',
        #p.adjust.method = 'BH',
        ref.group = 'uninfected')

  ## global FDR adjustment

  time_course$post_hoc_test <- time_course$post_hoc_test %>%
    compress(name = 'metabolite') %>%
    mutate(p.adj = p.adjust(p, 'BH'),
           metabolite = factor(metabolite,
                               names(time_course$post_hoc_test))) %>%
    blast(metabolite)

# Testing by Kruskal-Wallis test --------

  insert_msg('Testing for differences between the timepoints')

  time_course$test <- time_course$analysis_tbl %>%
    compare_variables(variables = globals$incov_lexicon$variable,
                      split_factor = 'timepoint',
                      what = 'test',
                      types = 'kruskal_test',
                      ci = FALSE,
                      exact = FALSE,
                      pub_styled = TRUE,
                      adj_method = 'BH')

# Plotting --------

  insert_msg('Plotting')

  time_course$plots <-
    list(variable = time_course$test$variable,
         plot_title = exchange(time_course$test$variable,
                               dict = globals$incov_lexicon,
                               key = 'variable',
                               value = 'label') %>%
           paste('INCOV', sep = ', '),
         plot_subtitle = paste('Kruskal-Wallis test',
                               time_course$test$significance)) %>%
    pmap(plot_variable,
         time_course$analysis_tbl,
         split_factor = 'timepoint',
         type = 'box',
         point_hjitter = 0,
         x_lab = 'Time after infection',
         y_lab = expression('normalized log'[2] * ' concentration'),
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_fill_manual(values = c('steelblue3',
                                       rep('indianred3', 6))) +
          theme(axis.title.x = element_blank())) %>%
    set_names(time_course$test$variable)

# Appending the plots with the post-hoc testing results ------

  insert_msg('Appending the plots with the post-hoc testing results')

  time_course$plots <- map2(time_course$plots,
                            time_course$post_hoc_test,
                            add_p_one)

# END ------

  insert_tail()
