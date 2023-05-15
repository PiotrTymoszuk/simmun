# Modeling of the metabolite responses by backward elimination

  insert_head()

# container ----

  mod_eli <- list()

# Construction of the null and full models -------

  insert_msg('construction of the null and full models')

  mod_eli$null_models <- list(formula = mod_globals$null_formulas) %>%
    pmap(make_lm,
         data = mod_globals$analysis_tbl,
         mod_fun = lm,
         family = NULL)

  mod_eli$full_models <-
    list(formula = mod_globals$full_formulas) %>%
    pmap(make_lm,
         data = mod_globals$analysis_tbl,
         mod_fun = lm,
         family = NULL)

# Backward elimination via step() -----

  insert_msg('Backward elimination')

  mod_eli$step_models <- mod_eli$full_models %>%
    map(step,
        step_fun = MASS::stepAIC,
        direction = 'backward',
        k = floor(log(nrow(stigma$analysis_tbl))))

# Checking assumptions of the optimized models -----

  insert_msg('Assumptions of the optimized models')

  ## normality and EOV of the residuals

  mod_eli$assumptions <- mod_eli$step_models %>%
    map(summary,
        type = 'assumptions',
        type.predict = 'response')

  ## diagnostic plots of the residuals

  mod_eli$resid_plots <- mod_eli$step_models %>%
    map(plot,
        type = 'residuals',
        cust_theme = globals$common_theme,
        type.predict = 'response')

  ## response - explanatory variable relationship

  mod_eli$relation_plots <- mod_eli$step_models %>%
    map(plot,
        type = 'relationship',
        cust_theme = globals$common_theme,
        type.predict = 'response')

# LRT versus the null model ------

  insert_msg('LRT versus the null model')

  mod_eli$lrt <- map2(mod_eli$step_models,
                        mod_eli$null_models,
                        ~anova(.x$model, .y$model))

  ## LRT summaries

  mod_eli$lrt_summary <- mod_eli$lrt %>%
    map(as.data.frame) %>%
    map(~.x[2, ]) %>%
    compress(names_to = 'response') %>%
    mutate(significance = ifelse(`Pr(>F)` < 0.05,
                                 paste('p =', signif(`Pr(>F)`, 2)),
                                 paste0('ns (p = ', signif(`Pr(>F)`, 2), ')')),
           plot_cap = paste0('F(', abs(Df), ', ',
                             Res.Df, ') = ', signif(`F`, 2)),
           plot_cap = paste(plot_cap, significance, sep = ', '),
           plot_cap = paste(plot_cap, Res.Df + 1, sep = ', n = '),
           plot_cap = set_names(plot_cap, response)) %>%
    as_tibble

# Fit goodness and errors ------

  insert_msg('Fit goodness and errors')

  mod_eli$fit_stats <- mod_eli$step_models %>%
    map_dfr(summary,
            type = 'fit')

# Cross-validation with caret ------

  insert_msg('Cross-validation')

  mod_eli$step_formulas <- mod_eli$step_models %>%
    map(formula)

  set.seed(123456)

  registerDoParallel(cores = 7)

  ## safely: there's no valid model for PHE

  mod_eli$caret_models <- mod_eli$step_formulas %>%
    map(~safely(train)(form = .x,
                       data = as.data.frame(mod_globals$analysis_tbl),
                       method = 'lm',
                       metric = 'RMSE',
                       trControl = trainControl(method = 'repeatedcv',
                                                number = 10,
                                                savePredictions = 'final',
                                                returnData = TRUE,
                                                returnResamp = 'final'))) %>%
    map(~.x$result) %>%
    compact %>%
    map(as_caretx)

  ## patching the training data, caret has problems with factors in formulas
  ## of linear models

  for(i in names(mod_eli$caret_models)) {

    mod_eli$caret_models[[i]]$trainingData <-
      mod_eli$caret_models[[i]]$trainingData %>%
      model.matrix(~., data = .) %>%
      as.data.frame

  }

  stopImplicitCluster()

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  plan('multisession')

  mod_eli$cv_stats <- mod_eli$caret_models %>%
    extract_caret_stats

  plan('sequential')

# Plots with the R-squares and RMSE for train and CV data ------

  insert_msg('Plots of the train and CV fit stats')

  mod_eli$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(mod_eli$cv_stats,
                     subtitle = 'Backward elimination') %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(mod_eli$caret_models)),
                           labels = exchange(names(mod_eli$caret_models),
                                             dict = globals$stigma_lexicon,
                                             key = 'variable',
                                             value = 'label')))

# Inference and Forest plots of the model betas -------

  insert_msg('Model inference and Forest plots of the betas')

  ## model estimates
  ## no valid model present for log_phe (intercept only)
  ## hence working with safely()

  mod_eli$inference <- mod_eli$step_models %>%
    map(safely(summary), type = 'inference') %>%
    map(~.x$result) %>%
    compact %>%
    map(mutate,
        var_lab = exchange(variable,
                           dict = globals$stigma_lexicon),
        var_lab = stri_replace(var_lab, fixed = ' score', replacement = ''),
        level = ifelse(level == 'CoV recovery' | level == 'yes',
                       '', level))

  mod_eli$forest_plots <-
    list(x = mod_eli$inference %>%
           map(filter, p_value < 0.05),
         plot_title = exchange(stigma$responses[names(mod_eli$inference)],
                               dict = globals$stigma_lexicon,
                               key = 'variable',
                               value = 'base_label') %>%
           paste('SIMMUN', sep = ', '),
         x_lab = exchange(stigma$responses[names(mod_eli$inference)],
                          dict = globals$stigma_lexicon) %>%
           paste0(', Z-score'),
         plot_subtitle = mod_eli$lrt_summary$plot_cap[names(mod_eli$inference)]) %>%
    pmap(plot_forest,
         variable = 'var_lab',
         hide_baseline = TRUE,
         estimate_size = 2.5,
         estimate_hjust = 0.5,
         cust_theme = globals$common_theme) %>%
    map(~.x +
          scale_y_discrete(labels = function(x) stri_replace(x, regex = ':\\s+\\n', replacement = '\n')))

# Appending the models with the baseline values ------

  insert_msg('Baseline values')

  mod_eli$intercepts <- mod_eli$inference[names(mod_eli$forest_plots)] %>%
    map_dfr(filter,
            level == 'baseline') %>%
    mutate(baseline_lab = paste0('baseline: ', signif(estimate, 2))) %>%
    .$baseline_lab

  mod_eli$forest_plots <-
    map2(mod_eli$forest_plots,
         mod_eli$intercepts,
         ~.x +
           labs(x = paste(.x$labels$x,
                          .y,
                          sep = ', ')))

# END ------

  rm(i)

  insert_tail()
