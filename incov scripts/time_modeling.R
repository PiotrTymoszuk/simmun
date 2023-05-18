# Additional analysis with robust linear modeling as requested by a reviewer

  insert_head()

# container ----

  time_rlm <- list()

# analysis globals ------

  insert_msg('Globals')

  ## variable lexicon

  time_rlm$lexicon <- globals$incov_lexicon %>%
    mutate(label_long = variable,
           variable = make.names(variable))

  ## analysis tables

  time_rlm$analysis_tbl <- incov$analysis_tbl %>%
    select(patient_id,
           timepoint,
           any_of(set_names(time_rlm$lexicon$label_long,
                            time_rlm$lexicon$variable))) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'uninfected'"))

  ## setting the baselines

  time_rlm$analysis_tbl <-
    list(uninfected = c('uninfected', 'acute', 'sub-acute', 'recovery'),
         acute = c('acute', 'uninfected', 'sub-acute', 'recovery')) %>%
    map(~mutate(time_rlm$analysis_tbl,
                timepoint = factor(timepoint, levels = .x)))

  ## model formulas and null model formulas

  time_rlm$formulas <- time_rlm$lexicon$variable %>%
    map(~paste0('`', .x, '` ~ timepoint')) %>%
    map(as.formula) %>%
    set_names(time_rlm$lexicon$variable)

  time_rlm$null_formulas <- time_rlm$lexicon$variable %>%
    map(~paste0('`', .x, '` ~ 1')) %>%
    map(as.formula) %>%
    set_names(time_rlm$lexicon$variable)

  ## colors for plots

  time_rlm$plot_colors <-
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

# Construction of the models -----

  insert_msg('Model construction')

  ## models and null models

  time_rlm$models <- time_rlm$analysis_tbl %>%
    map(function(data) time_rlm$formulas %>%
          map(~make_lm(data = data,
                       formula = .x,
                       mod_fun = rlm,
                       family = NULL,
                       method = 'MM',
                       psi = psi.huber)))

# Inference --------

  insert_msg('Inference')

  ## and global FDR correction

  time_rlm$inference <- time_rlm$models %>%
    map(~map(.x,
             summary,
             type = 'inference',
             ci_method = 'distribution')) %>%
    map(compress,
        names_to = 'metabolite') %>%
    map(re_adjust) %>%
    map(mutate,
        metabolite = factor(metabolite, names(time_rlm$models[[1]]))) %>%
    map(blast, metabolite)

# Plotting ------

  insert_msg('Plotting')

  time_rlm$plots <-
    list(variable = time_rlm$lexicon$variable,
         plot_title = exchange(time_rlm$lexicon$variable,
                               dict = time_rlm$lexicon,
                               key = 'variable',
                               value = 'label') %>%
           paste('INCOV', sep = ', ')) %>%
    pmap(plot_variable,
         time_rlm$analysis_tbl$uninfected,
         split_factor = 'timepoint',
         type = 'box',
         point_hjitter = 0,
         x_lab = 'Time after infection',
         y_lab = expression('normalized log'[2] * ' concentration'),
         cust_theme = globals$common_theme,
         x_n_labs = TRUE) %>%
    map(~.x +
          scale_fill_manual(values = c('steelblue3',
                                       rep('indianred3', 6))) +
          theme(axis.title.x = element_blank())) %>%
    set_names(time_course$test$variable)

# Appending the plots with model betas and p values ------

  insert_msg('Appending the plots with model betas and p values')

  time_rlm$plots <-
    map2(time_rlm$plots,
         time_rlm$inference$uninfected,
         add_beta)

# Ribbon plots for the inflammatory and metabolic variables -----

  insert_msg('Ribbon plots')

  ## ribbon plots for inflammatory cytokines, TRP-metabolites and PHE/TYR

  time_rlm$ribbon_plots <-
    list(vars = list(inflammatory = c('IL6_INF', 'IL10_INF', 'TNF_INF', 'IFNG_INF'),
                     trp = c('tryptophan', 'kynurenine', 'quinolinate', 'serotonin'),
                     tyr = c('phenylalanine', 'tyrosine', 'dopamine.3.O.sulfate')),
         plot_title = paste0(c('Inflammatory cytokines',
                               'TRP decay and 5-HT',
                               'PHE/TYR and DA'),
                             ', INCOV')) %>%
    pmap(plot_tc_ribbon,
         data = time_rlm$analysis_tbl$uninfected,
         plot_subtitle = 'Mean, 2 \u00D7 SEM',
         average_stat = 'mean',
         ribbon_stat = '2SEM',
         ribbon_alpha = 0.3) %>%
    map2(., time_rlm$plot_colors ,
         ~.x +
           scale_color_manual(values = .y,
                              labels = function(x) exchange(x, dict = time_rlm$lexicon)) +
           scale_fill_manual(values = .y,
                             labels = function(x) exchange(x, dict = time_rlm$lexicon)))

# Forest plots with model betas and their 95% CI ------

  insert_msg('Forest plots')

  time_rlm$forest_plots <-
    list(x = time_rlm$inference,
         y = paste('Robust linear modeling,',
                   c('uninfected', 'acute')),
         z = c('uninfected', 'acute')) %>%
    pmap(function(x, y, z) list(vars = list(inflammatory = c('IL6_INF',
                                                             'IL10_INF',
                                                             'TNF_INF',
                                                             'IFNG_INF'),
                                            trp = c('tryptophan',
                                                    'kynurenine',
                                                    'quinolinate',
                                                    'serotonin'),
                                            tyr = c('phenylalanine',
                                                    'tyrosine',
                                                    'dopamine.3.O.sulfate')),
                                plot_title = paste(c('Inflammatory cytokines',
                                                     'TRP decay and 5-HT',
                                                     'PHE/TYR and DA'))) %>%
           pmap(plot_tc_beta,
                data = x,
                plot_subtitle = y,
                baseline_lab = z,
                hide_baseline = FALSE,
                show_connector = FALSE) %>%
           map2(., time_rlm$plot_colors,
                ~.x +
                  scale_color_manual(values = .y,
                                     labels = function(x) exchange(x,
                                                                   dict = time_rlm$lexicon)) +
                  scale_fill_manual(values = .y,
                                    labels = function(x) exchange(x,
                                                                  dict = time_rlm$lexicon))))


# Tables with the modeling results ------

  insert_msg('Tables with the inference results')

  time_rlm$result_tbl <- time_rlm$inference %>%
    map(~map(.x,
             ~mutate(.x,
                     n = ifelse(!is.na(n),
                                n = .x$n_complete[[1]] - sum(n, na.rm = TRUE),
                                n))) %>%
          compress(names_to = 'response') %>%
          transmute(Response = exchange(response, dict = time_rlm$lexicon),
                    Timepoint = level,
                    `Estimate, 95% CI` = paste0(signif(estimate, 2), ' [',
                                                signif(lower_ci, 2), ' - ',
                                                signif(upper_ci, 2), ']'),
                    Significance = significance)) %>%
    map2(., c('uninfected', 'acute'),
         ~mutate(.x, Baseline = .y)) %>%
    map(relocate, Baseline)

# END ------

  insert_tail()
