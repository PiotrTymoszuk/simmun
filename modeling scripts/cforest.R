# Modeling metabolite responses with ranger's Random forest

  insert_head()

# container ----

  mod_crf <- list()

# analysis globals -----

  insert_msg('Analysis globals')

  ## since there's no normality assumption, variables
  ## are fed into the model without normality-stabilizing transformations

  mod_crf$lexicon <- stigma$expl_lexicon %>%
    mutate(variable = stri_replace(variable,
                                   regex = '^(log|sqrt)_',
                                   replacement = ''))

  mod_crf$responses <- stigma$responses %>%
    stri_replace(regex = '^(log|sqrt)_',
                 replacement = '') %>%
    set_names(names(stigma$responses))

  mod_crf$analysis_tbl <- stigma$data %>%
    mutate(infection = car::recode(cov,
                                   "'healthy' = 'uninfected';
                             'SARS-CoV-2' = 'recovery'"),
           infection = factor(infection,
                              c('uninfected', 'recovery')),
           age = age/10) %>%
    select(patient_id,
           all_of(mod_crf$responses),
           all_of(mod_crf$lexicon$variable)) %>%
    filter(complete.cases(.))

  ## normalizing the explanatory variables

  mod_crf$analysis_tbl[mod_crf$lexicon$variable] <-
    mod_crf$analysis_tbl[mod_crf$lexicon$variable] %>%
    map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x)

  mod_crf$analysis_tbl <- mod_crf$analysis_tbl %>%
    column_to_rownames('patient_id')

  ## model formulas

  mod_crf$formulas <- mod_crf$responses %>%
    map(~paste(.x,
               paste(mod_crf$lexicon$variable, collapse = ' + '),
               sep = ' ~ ')) %>%
    map(as.formula)


  ## tuning grids

  mod_crf$tune_grid <-
    expand.grid(mtry = c(2, 4, 6, 12))

# Caret models -----

  insert_msg('Caret models')

  set.seed(123456)

  registerDoParallel(cores = 7)

  mod_crf$caret_models <- mod_crf$formulas %>%
    map(~train(form = .x,
               data = mod_crf$analysis_tbl,
               method = 'cforest',
               metric = 'RMSE',
               tuneGrid =  mod_crf$tune_grid,
               trControl = trainControl(method = 'repeatedcv',
                                        number = 10,
                                        savePredictions = 'final',
                                        returnData = TRUE,
                                        returnResamp = 'final'),
               controls = cforest_unbiased(ntree = 1000))) %>%
    map(as_caretx)

  stopImplicitCluster()

# Futures -----

  insert_msg('Registering the futures')

  plan('multisession')

# Plots of model residuals -----

  insert_msg('Plots of model residuals')

  mod_crf$resid_plots <- mod_crf$caret_models %>%
    future_map(plot,
               cust_theme = globals$common_theme,
               .options = furrr_options(seed = TRUE))

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  mod_crf$cv_stats <- mod_crf$caret_models %>%
    extract_caret_stats

# Plotting the R-squares and RMSE for the training and CV datasets ------

  insert_msg('Plotting the data and CV stats')

  mod_crf$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(mod_crf$cv_stats,
                     subtitle = 'Conditional Random Forest') %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(mod_crf$caret_models)),
                     labels = exchange(names(mod_crf$caret_models),
                                       dict = globals$stigma_lexicon,
                                       key = 'variable',
                                       value = 'label')))

# Variable importance -----

  insert_msg('Variable importance')

  ## extraction regex

  mod_crf$importance$ext_regex <- mod_crf$lexicon$variable %>%
    sort(decreasing = TRUE) %>%
    paste(collapse = '|')

  ## permutation variable importance

  mod_crf$importance$stats <- mod_crf$caret_models %>%
    map(~.x$finalModel) %>%
    map(varimp) %>%
    map(compress,
        names_to = 'parameter',
        values_to = 'importance') %>%
    map(mutate,
        variable = stri_extract(parameter,
                                regex = mod_crf$importance$ext_regex),
        level = stri_replace(parameter,
                             regex = mod_crf$importance$ext_regex,
                             replacement = ''),
        axis_lab = exchange(variable,
                            dict = mod_crf$lexicon),
        axis_lab = ifelse(level == '' | level == 'yes',
                          axis_lab,
                          paste(axis_lab, level, sep = ': ')))

  ## plots

  mod_crf$importance$plots <- mod_crf$importance$stats %>%
    map2(.,
         exchange(names(.), dict = globals$stigma_lexicon),
         ~.x %>%
           ggplot(aes(x = importance,
                      y = reorder(axis_lab, importance))) +
           geom_bar(color = 'black',
                    fill = 'steelblue',
                    stat = 'identity') +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = .y,
                subtitle = 'Permutation importance, Conditional Random Forest',
                x = 'importance'))

# END -----

  plan('sequential')

  insert_tail()
