# Testing of effects of BMI on the time course of metabolites and cytokines
# Analyzed with robust linear modeling, I'm skipping the interaction effect
# since it was not significant for any of the responses

  insert_head()

# container ----

  bmi_incov <- list()

# analysis globals ------

  insert_msg('Analysis globals')

  ## variable lexicon

  bmi_incov$lexicon <- globals$incov_lexicon %>%
    mutate(label_long = variable,
           variable = make.names(variable))

  ## analysis tables

  bmi_incov$analysis_tbl <- incov$analysis_tbl %>%
    select(patient_id,
           timepoint,
           bmi_class,
           any_of(set_names(bmi_incov$lexicon$label_long,
                            bmi_incov$lexicon$variable))) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"),
           timepoint = factor(timepoint,
                              c('uninfected', 'acute', 'sub-acute', 'recovery')))

  ## colors for plots

  bmi_incov$plot_colors <-
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

  bmi_incov$plot_titles <-
    paste0(c('Inflammatory cytokines',
             'Serotonin and kynurenine pathway',
             'Catecholamine pathway'),
           ', INCOV') %>%
    set_names(c('inflammatory', 'trp', 'tyr'))

  ## model formulas

  bmi_incov$formulas <- bmi_incov$lexicon$variable %>%
    map(~paste(.x, '~ bmi_class + timepoint')) %>%
    map(as.formula) %>%
    set_names(bmi_incov$lexicon$variable)

# Construction of the models ------

  insert_msg('Construction of the models')

  bmi_incov$models <- bmi_incov$formulas %>%
    map(~make_lm(data = bmi_incov$analysis_tbl,
                 formula = .x,
                 mod_fun = rlm,
                 family = NULL,
                 method = 'MM',
                 psi = psi.huber))

# Model inference -----

  insert_msg('Model inference')

  bmi_incov$inference <- bmi_incov$models %>%
    map(summary, 'inference', ci_method = 'distribution') %>%
    map(re_adjust, method = 'none') %>%
    map(mutate,
        variable = stri_replace(variable,
                                fixed = 'bmi_class',
                                replacement = 'BMI'),
        plot_lab = paste0(signif(estimate, 2), ' [',
                          signif(lower_ci, 2), ' - ',
                          signif(upper_ci, 2), ']'),
        var_lab = paste(variable, level, sep = ': '),
        plot_lab = paste(var_lab, plot_lab, sep = ': '),
        plot_lab = paste(plot_lab, significance, sep = ', '))

  bmi_incov$inference_tags <- bmi_incov$inference %>%
    map(~.x$plot_lab[2:5]) %>%
    map(paste, collapse = '\n')

# Ribbon plots, spit by timepoint -------

  insert_msg('Ribbon plots, spit be the response type')

  for(i in levels(bmi_incov$analysis_tbl$bmi_class)) {

    bmi_incov$ribbon_plots[[i]] <-
      list(vars = list(inflammatory = c('IL6_INF', 'IL10_INF', 'TNF_INF', 'IFNG_INF'),
                       trp = c('tryptophan', 'kynurenine', 'quinolinate', 'serotonin'),
                       tyr = c('phenylalanine', 'tyrosine', 'dopamine.3.O.sulfate')),
           plot_title = paste(bmi_incov$plot_titles, i, sep = ', ')) %>%
      pmap(plot_tc_ribbon,
           data = blast(bmi_incov$analysis_tbl, bmi_class)[[i]],
           plot_subtitle = 'Mean, 2 \u00D7 SEM',
           average_stat = 'mean',
           ribbon_stat = '2SEM',
           ribbon_alpha = 0.3) %>%
      map2(., bmi_incov$plot_colors ,
           ~.x +
             scale_color_manual(values = .y,
                                labels = function(x) exchange(x, dict = bmi_incov$lexicon)) +
             scale_fill_manual(values = .y,
                               labels = function(x) exchange(x, dict = bmi_incov$lexicon)))


  }

# Box plots split by the BMI class -----

  insert_msg('Box plots, split by the BMI class')

  bmi_incov$plots <-
    list(x = bmi_incov$lexicon$variable,
         y = bmi_incov$lexicon$label,
         z = bmi_incov$inference_tags) %>%
    pmap(function(x, y, z) bmi_incov$analysis_tbl %>%
           ggplot(aes(x = timepoint,
                      y = .data[[x]],
                      fill = bmi_class)) +
           geom_boxplot(outlier.colour = NA,
                        alpha = 0.25,
                        position = position_dodge(0.8)) +
           geom_point(shape = 21,
                      size = 2,
                      position = position_jitterdodge(jitter.width = 0.1,
                                                      jitter.height = 0,
                                                      dodge.width = 0.8)) +
           scale_fill_manual(values = c(normal = 'steelblue',
                                        overweight = 'coral3',
                                        obesity = 'coral4'),
                             name = '') +
           globals$common_theme +
           labs(title = y,
                y = expression('normalized log'[2] * ' concentration'),
                x = 'SARS-CoV-2 infection',
                tag = paste0('\n', z))) %>%
    set_names(bmi_incov$lexicon$variable)

# END -----

  rm(i)

  insert_tail()
