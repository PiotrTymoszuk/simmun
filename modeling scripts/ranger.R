# Modeling metabolite responses with ranger's Random forest

  insert_head()

# container ----

  mod_ranger <- list()

# analysis globals -----

  insert_msg('Analysis globals')

  ## since there's no normality assumption, variables
  ## are fed into the model without normality-stabilizing transformations

  mod_ranger$lexicon <- stigma$expl_lexicon %>%
    mutate(variable = stri_replace(variable,
                                   regex = '^(log|sqrt)_',
                                   replacement = ''))

  mod_ranger$responses <- stigma$responses %>%
    stri_replace(regex = '^(log|sqrt)_',
                 replacement = '') %>%
    set_names(names(stigma$responses))

  mod_ranger$analysis_tbl <- stigma$data %>%
    mutate(infection = car::recode(cov,
                                   "'healthy' = 'uninfected';
                             'SARS-CoV-2' = 'recovery'"),
           infection = factor(infection,
                              c('uninfected', 'recovery')),
           age = age/10) %>%
    select(patient_id,
           all_of(mod_ranger$responses),
           all_of(mod_ranger$lexicon$variable)) %>%
    filter(complete.cases(.))

  ## normalizing the explanatory variables

  mod_ranger$analysis_tbl[mod_ranger$lexicon$variable] <-
    mod_ranger$analysis_tbl[mod_ranger$lexicon$variable] %>%
    map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x)

  mod_ranger$analysis_tbl <- mod_ranger$analysis_tbl %>%
    column_to_rownames('patient_id')

  ## model formulas

  mod_ranger$formulas <- mod_ranger$responses %>%
    map(~paste(.x,
               paste(mod_ranger$lexicon$variable, collapse = ' + '),
               sep = ' ~ ')) %>%
    map(as.formula)


  ## tuning grids

  mod_ranger$tune_grid <-
    expand.grid(mtry = c(2, 4, 6, 12),
                splitrule = c('variance', 'extratrees', 'maxstat'),
                min.node.size = c(3, 5))

# Caret models -----

  insert_msg('Caret models')

  set.seed(123456)

  registerDoParallel(cores = 7)

  mod_ranger$caret_models <- mod_ranger$formulas %>%
    map(~train(form = .x,
               data = mod_ranger$analysis_tbl,
               method = 'ranger',
               metric = 'RMSE',
               tuneGrid =  mod_ranger$tune_grid,
               trControl = trainControl(method = 'repeatedcv',
                                        number = 10,
                                        savePredictions = 'final',
                                        returnData = TRUE,
                                        returnResamp = 'final')),
        num.trees = 1000) %>%
    map(as_caretx)

  stopImplicitCluster()

# Futures -----

  insert_msg('Registering the futures')

  plan('multisession')

# Plots of model residuals -----

  insert_msg('Plots of model residuals')

  mod_ranger$resid_plots <- mod_ranger$caret_models %>%
    future_map(plot,
               cust_theme = globals$common_theme,
               .options = furrr_options(seed = TRUE))

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  mod_ranger$cv_stats <- mod_ranger$caret_models %>%
    extract_caret_stats

# Plotting the R-squares and RMSE for the training and CV datasets ------

  insert_msg('Plotting the data and CV stats')

  mod_ranger$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(mod_ranger$cv_stats,
                     subtitle = 'Random Forest') %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(mod_ranger$caret_models)),
                           labels = exchange(names(mod_ranger$caret_models),
                                             dict = globals$stigma_lexicon,
                                             key = 'variable',
                                             value = 'label')))

# Variable importance -----

  insert_msg('Variable importance')

  mod_ranger$importance$models <-
    list(x = mod_ranger$formulas,
         y = mod_ranger$caret_models) %>%
    pmap(function(x, y) ranger(formula = x,
                               data = mod_ranger$analysis_tbl,
                               importance = 'permutation',
                               num.trees = 1000,
                               mtry = y$bestTune$mtry,
                               splitrule = y$bestTune$splitrule,
                               min.node.size = y$bestTune$min.mod.size))

  mod_ranger$importance$stats <- mod_ranger$importance$models %>%
    map(importance) %>%
    map(compress,
        names_to = 'variable',
        values_to = 'importance')

  mod_ranger$importance$plots <- mod_ranger$importance$stats %>%
    map2(.,
         exchange(names(.), dict = globals$stigma_lexicon),
         ~.x %>%
           ggplot(aes(x = importance,
                      y = reorder(variable, importance))) +
           geom_bar(color = 'black',
                    fill = 'steelblue',
                    stat = 'identity') +
           scale_y_discrete(labels = function(x) exchange(x, dict = globals$stigma_lexicon)) +
           globals$common_theme +
           theme(axis.title.y = element_blank()) +
           labs(title = .y,
                subtitle = 'Permutation importance, Random Forest',
                x = 'importance'))


# END -----

  plan('sequential')

  insert_tail()
