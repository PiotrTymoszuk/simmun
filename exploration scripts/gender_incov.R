# Testing of effects of gender on the time course of metabolites and cytokines
# Analyzed with robust linear modeling, I'm skipping the interaction effect
# since it was not significant for any of the responses


  insert_head()

# container ----

  gender_incov <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variable lexicon

  gender_incov$lexicon <- globals$incov_lexicon %>%
    mutate(label_long = variable,
           variable = make.names(variable))

  ## analysis tables

  gender_incov$analysis_tbl <- incov$analysis_tbl %>%
    select(patient_id,
           timepoint,
           sex,
           any_of(set_names(gender_incov$lexicon$label_long,
                            gender_incov$lexicon$variable))) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"),
           timepoint = factor(timepoint,
                              c('uninfected', 'acute', 'sub-acute', 'recovery')))

  ## colors for plots

  gender_incov$plot_colors <-
    list(inflammatory = c('IL6_INF' = 'orangered2',
                          'IL10_INF' = 'orangered4',
                          'TNF_INF' = 'gray60',
                          'IFNG_INF' = 'black'),
         trp = c('tryptophan' = 'coral',
                 'kynurenine' = 'indianred2',
                 'quinolinate' = 'indianred4',
                 'serotonin' = 'steelblue'),
         tyr = c('phenylalanine' = 'darkolivegreen3',
                 'tyrosine' = 'darkolivegreen4',
                 'dopamine.3.O.sulfate' = 'steelblue'))

  ## plot titles

  gender_incov$plot_titles <-
    paste0(c('Inflammatory cytokines',
             'Serotonin and kynurenine pathway',
             'Catecholamine pathway'),
           ', INCOV') %>%
    set_names(c('inflammatory', 'trp', 'tyr'))

  ## model formulas

  gender_incov$formulas <- gender_incov$lexicon$variable %>%
    map(~paste(.x, '~ sex + timepoint')) %>%
    map(as.formula) %>%
    set_names(gender_incov$lexicon$variable)

# Construction of the models ------

  insert_msg('Construction of the models')

  gender_incov$models <- gender_incov$formulas %>%
    map(~make_lm(data = gender_incov$analysis_tbl,
                 formula = .x,
                 mod_fun = rlm,
                 family = NULL,
                 method = 'MM',
                 psi = psi.huber))

# Model inference -----

  insert_msg('Model inference')

  gender_incov$inference <- gender_incov$models %>%
    map(summary, 'inference', ci_method = 'distribution') %>%
    map(re_adjust, method = 'none') %>%
    map(mutate,
        plot_lab = paste0(signif(estimate, 2), ' [',
                          signif(lower_ci, 2), ' - ',
                          signif(upper_ci, 2), ']'),
        var_lab = paste(variable, level, sep = ': '),
        plot_lab = paste(var_lab, plot_lab, sep = ': '),
        plot_lab = paste(plot_lab, significance, sep = ', '))

  gender_incov$inference_tags <- gender_incov$inference %>%
    map(~.x$plot_lab[2:5]) %>%
    map(paste, collapse = '\n')

# Ribbon plots, spit by timepoint -------

  insert_msg('Ribbon plots, spit be the response type')

  for(i in levels(gender_incov$analysis_tbl$sex)) {

    gender_incov$ribbon_plots[[i]] <-
      list(vars = list(inflammatory = c('IL6_INF', 'IL10_INF', 'TNF_INF', 'IFNG_INF'),
                       trp = c('tryptophan', 'kynurenine', 'quinolinate', 'serotonin'),
                       tyr = c('phenylalanine', 'tyrosine', 'dopamine.3.O.sulfate')),
           plot_title = paste(gender_incov$plot_titles, i, sep = ', ')) %>%
      pmap(plot_tc_ribbon,
           data = blast(gender_incov$analysis_tbl, sex)[[i]],
           plot_subtitle = 'Mean, 2 \u00D7 SEM',
           average_stat = 'mean',
           ribbon_stat = '2SEM',
           ribbon_alpha = 0.3) %>%
      map2(., gender_incov$plot_colors,
           ~.x +
             scale_color_manual(values = .y,
                                labels = function(x) exchange(x, dict = time_rlm$lexicon)) +
             scale_fill_manual(values = .y,
                               labels = function(x) exchange(x, dict = time_rlm$lexicon)))


  }

  gender_incov$ribbon_plots <- gender_incov$ribbon_plots %>%
    transpose %>%
    map(~set_common_y_axis_range(.x[[1]], .x[[2]])) %>%
    map(set_names, levels(gender_incov$analysis_tbl$sex))

# Box plots split by the gender -----

  insert_msg('Box plots, split by the gender')

  gender_incov$plots <-
    list(x = gender_incov$lexicon$variable,
         y = gender_incov$lexicon$label,
         z = gender_incov$inference_tags) %>%
    pmap(function(x, y, z) gender_incov$analysis_tbl %>%
           ggplot(aes(x = timepoint,
                      y = .data[[x]],
                      fill = sex)) +
           geom_boxplot(outlier.colour = NA,
                        alpha = 0.25,
                        position = position_dodge(0.8)) +
           geom_point(shape = 21,
                      size = 2,
                      position = position_jitterdodge(jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      dodge.width = 0.8)) +
           scale_fill_manual(values = c(female = 'steelblue',
                                        male = 'indianred3'),
                             name = '') +
           globals$common_theme +
           labs(title = y,
                y = expression('normalized log'[2] * ' concentration'),
                x = 'SARS-CoV-2 infection',
                tag = paste0('\n', z))) %>%
    set_names(gender_incov$lexicon$variable)

# END -----

  rm(i)

  insert_tail()
