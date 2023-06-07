# Comparison of metabolite and inflammation markers in the BMI classes
# and CoV infection. Done with two-way ANOVA with the interaction term.

  insert_head()

# container ------

  bmi_simmun <- list()

# analysis globals -------

  insert_msg('Analysis globals')

  ## variables

  bmi_simmun$variables <- c(neo = 'log_neo',
                               nlr = 'log_nlr',
                               stigma$responses)

  ## analysis table

  bmi_simmun$analysis_tbl <-
    stigma$data[c('patient_id', 'bmi_class', 'cov', bmi_simmun$variables)] %>%
    mutate(cov = car::recode(cov, "'healthy' = 'uninfected'"),
           cov = factor(cov, c('uninfected', 'SARS-CoV-2'))) %>%
    filter(!is.na(bmi_class))

  ## model formulas

  bmi_simmun$formulas <- bmi_simmun$variables %>%
    map(~paste(.x, '~ bmi_class * cov')) %>%
    map(as.formula)

# Models -----

  insert_msg('Models')

  bmi_simmun$models <- bmi_simmun$formulas %>%
    map(~make_lm(data = bmi_simmun$analysis_tbl,
                 formula = .x,
                 mod_fun = lm,
                 family = NULL))

# Assumptions -----

  insert_msg('Assumptions')

  ## there moderate violations of the normality assumption for the
  ## model residuals measured by Shapiro-Wilk test and QQ plots of the residuals
  ## (neopterin, kyn/trp, tyr, phe/tyr)

  bmi_simmun$assumptions <- bmi_simmun$models %>%
    map(summary, 'assumptions', type.predict = 'response')

  bmi_simmun$resid_plots <- bmi_simmun$models %>%
    map(plot, cust_theme = globals$common_theme, type.predict = 'response')

# Fit stats ------

  insert_msg('Fit stats')

  bmi_simmun$fit_stats <- bmi_simmun$models %>%
    map(summary, 'fit') %>%
    map(mutate, transformation = 'response') %>%
    compress(names_to = 'response')

# ANOVA ------

  insert_msg('Two-way ANOVA')

  bmi_simmun$test <- bmi_simmun$models %>%
    map(anova) %>%
    map(re_adjust, p_variable = 'Pr(>F)', method = 'none') %>%
    map(~mutate(.x,
                variable = stri_replace(variable,
                                        fixed = 'cov',
                                        replacement = 'SARS-CoV-2'),
                variable = stri_replace(variable,
                                        fixed = 'bmi_class',
                                        replacement = 'BMI'),
                df_sum = paste(Df, .x$Df[4], sep = ', '),
                eta_lab = paste('\u03B7\u00B2 =',
                                signif(frac_explained, 2)),
                eta_lab = paste(variable, eta_lab, significance, sep = ', ')))

  bmi_simmun$plot_tags <- bmi_simmun$test %>%
    map(~.x$eta_lab[1:3]) %>%
    map_chr(paste, collapse = '\n')

# Plotting ------

  insert_msg('Plotting')

  bmi_simmun$plots <-
    list(x = bmi_simmun$variables,
         y = exchange(names(bmi_simmun$variables),
                      globals$stigma_lexicon),
         z = exchange(bmi_simmun$variables,
                      globals$stigma_lexicon,
                      value = 'axis_label'),
         v = bmi_simmun$plot_tags) %>%
    pmap(function(x, y, z, v) bmi_simmun$analysis_tbl %>%
           ggplot(aes(x = cov,
                      y = .data[[x]],
                      fill = bmi_class)) +
           geom_boxplot(outlier.color = NA,
                        alpha = 0.25,
                        position = position_dodge(0.9)) +
           geom_point(shape = 21,
                      size = 2,
                      alpha = 0.8,
                      position = position_jitterdodge(jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      dodge.width = 0.9)) +
           scale_fill_manual(values = c(normal = 'steelblue',
                                        overweight = 'coral3',
                                        obesity = 'coral4'),
                             name = '') +
           globals$common_theme +
           theme(axis.title.x = element_blank()) +
           labs(title = y,
                y = z,
                tag = paste0('\n', v)))

# END -----

  insert_tail()
