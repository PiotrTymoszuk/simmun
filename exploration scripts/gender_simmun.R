# Comparison of metabolite and inflammation markers in the genders
# and CoV infection. Done with two-way ANOVA with the interaction term.

  insert_head()

# container ------

  gender_simmun <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## variables

  gender_simmun$variables <- c(neo = 'log_neo',
                               nlr = 'log_nlr',
                               stigma$responses)

  ## analysis table

  gender_simmun$analysis_tbl <-
    stigma$data[c('patient_id', 'sex', 'cov', gender_simmun$variables)] %>%
    mutate(cov = car::recode(cov, "'healthy' = 'uninfected'"),
           cov = factor(cov, c('uninfected', 'SARS-CoV-2')))

  ## model formulas

  gender_simmun$formulas <- gender_simmun$variables %>%
    map(~paste(.x, '~ sex * cov')) %>%
    map(as.formula)

# Models -----

  insert_msg('Models')

  gender_simmun$models <- gender_simmun$formulas %>%
    map(~make_lm(data = gender_simmun$analysis_tbl,
                 formula = .x,
                 mod_fun = lm,
                 family = NULL))

# Assumptions -----

  insert_msg('Assumptions')

  ## there moderate violations of the normality assumption for the
  ## model residuals measured by Shapiro-Wilk test and QQ plots of the residuals
  ## (neopterin, tyrosine and phe/tyr ratio)

  gender_simmun$assumptions <- gender_simmun$models %>%
    map(summary, 'assumptions', type.predict = 'response')

  gender_simmun$resid_plots <- gender_simmun$models %>%
    map(plot, cust_theme = globals$common_theme, type.predict = 'response')

# Fit stats ------

  insert_msg('Fit stats')

  gender_simmun$fit_stats <- gender_simmun$models %>%
    map(summary, 'fit') %>%
    map(mutate, transformation = 'response') %>%
    compress(names_to = 'response')

# ANOVA ------

  insert_msg('Two-way ANOVA')

  gender_simmun$test <- gender_simmun$models %>%
    map(anova) %>%
    map(re_adjust, p_variable = 'Pr(>F)', method = 'none') %>%
    map(~mutate(.x,
                variable = stri_replace(variable,
                                        fixed = 'cov',
                                        replacement = 'SARS-CoV-2'),
                df_sum = paste(Df, .x$Df[4], sep = ', '),
                eta_lab = paste('\u03B7\u00B2 =',
                                signif(frac_explained, 2)),
                eta_lab = paste(variable, eta_lab, significance, sep = ', ')))

  gender_simmun$plot_tags <- gender_simmun$test %>%
    map(~.x$eta_lab[1:3]) %>%
    map_chr(paste, collapse = '\n')

# Plotting ------

  insert_msg('Plotting')

  gender_simmun$plots <-
    list(x = gender_simmun$variables,
         y = exchange(names(gender_simmun$variables),
                      globals$stigma_lexicon),
         z = exchange(gender_simmun$variables,
                      globals$stigma_lexicon,
                      value = 'axis_label'),
         v = gender_simmun$plot_tags) %>%
    pmap(function(x, y, z, v) gender_simmun$analysis_tbl %>%
           ggplot(aes(x = cov,
                      y = .data[[x]],
                      fill = sex)) +
           geom_boxplot(outlier.color = NA,
                        alpha = 0.25,
                        position = position_dodge(0.9)) +
           geom_point(shape = 21,
                      size = 2,
                      alpha = 0.8,
                      position = position_jitterdodge(jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      dodge.width = 0.9)) +
           scale_fill_manual(values = c(male = 'indianred3',
                                        female = 'steelblue'),
                             name = '') +
           globals$common_theme +
           theme(axis.title.x = element_blank()) +
           labs(title = y,
                y = z,
                tag = paste0('\n', v)))

# END -----

  insert_tail()
