# Modeling of the metabolite responses with Elastic Net regression

  insert_head()

# container ------

  mod_net <- list()

# analysis globals -----

  insert_msg('Analysis globals')

  ## model matrices

  mod_net$y <- mod_globals$analysis_tbl[stigma$responses] %>%
    as.list %>%
    set_names(names(stigma$responses))

  mod_net$x <-
    mod_globals$analysis_tbl[c('patient_id', stigma$expl_lexicon$variable)] %>%
    column_to_rownames('patient_id') %>%
    model.matrix(~., data = .)

  ## model constants

  mod_net$alpha <- 0.5

  mod_net$repetitions <- 1:100 %>%
    set_names(paste0('rep_', 1:100))

  ## CV folds

  set.seed(1234)

  mod_net$cv_folds <- mod_net$repetitions %>%
    map(function(x) createFolds(y = stigma$analysis_tbl$patient_id,
                                k = 100,
                                list = FALSE,
                                returnTrain = TRUE))

# Tuning of the glmnet model -----

  insert_msg('Finding the optimal lambda')

  ## models

  plan('multisession')

  mod_net$tuning$models <- mod_net$y %>%
    map(function(resp) mod_net$cv_folds %>%
          future_map(~cv.glmnet(x = mod_net$x,
                                y = resp,
                                alpha = mod_net$alpha,
                                type.measure = 'deviance',
                                foldid = .x,
                                family = 'gaussian'),
                     .options = furrr_options(seed = TRUE)))

  plan('sequential')

  ## lambdas

  mod_net$tuning$stats <- mod_net$tuning$models %>%
    map(~map(.x,
             ~filter(as_tibble(.x[c('lambda', 'nzero', 'cvm', 'cvlo', 'cvup')]),
                     lambda == .x[['lambda.min']]))) %>%
    map(compress,
        names_to = 'repetition') %>%
    map(arrange, cvm) %>%
    map(mutate,
        optimum = ifelse(cvm == min(cvm), 'yes', 'no'),
        plot_lab = ifelse(optimum == 'yes',
                          signif(lambda, 3),
                          NA))

  ## optimum per response

  mod_net$tuning$opt_lambda <- mod_net$tuning$stats %>%
    map(filter, optimum == 'yes') %>%
    map(~.x$lambda[1])

  ## plotting

  mod_net$tuning$cvm_plot <- mod_net$tuning$stats %>%
    map2(., exchange(names(.), dict = globals$stigma_lexicon),
         ~.x %>%
           ggplot(aes(x = lambda,
                      y = cvm,
                      color = optimum)) +
           geom_point(size = 2,
                      shape = 16) +
           geom_text_repel(aes(label = plot_lab),
                           size = 2.75) +
           scale_color_manual(values = c(no = 'steelblue',
                                         yes = 'coral3'),
                              name = 'Optimum') +
           globals$common_theme +
           labs(title = .y,
                subtitle = 'LASSO lambda tuning',
                x = expression(lambda[LASSO]),
                y = 'Deviance'))

  mod_net$tuning$df_plot <- mod_net$tuning$stats %>%
    map2(., exchange(names(.), dict = globals$stigma_lexicon),
         ~.x %>%
           ggplot(aes(x = nzero,
                      y = lambda,
                      color = optimum)) +
           geom_point(size = 2,
                      shape = 16) +
           geom_text_repel(aes(label = plot_lab),
                           size = 2.75) +
           scale_color_manual(values = c(no = 'steelblue',
                                         yes = 'coral3'),
                              name = 'Optimum') +
           globals$common_theme +
           labs(title = .y,
                subtitle = 'LASSO lambda tuning',
                y = expression(lambda[LASSO]),
                x = 'Non-zero coefficients'))

# Construction of GLMNET models with the final lambda ------

  insert_msg('GLMNET models')

  mod_net$glmnet_models <-
    list(y = mod_net$y,
         lambda = mod_net$tuning$opt_lambda) %>%
    pmap(glmnet,
         x = mod_net$x,
         family = 'gaussian',
         alpha = mod_net$alpha)

# Cross-validation with caret ------

  insert_msg('Cross-validation')

  set.seed(123456)

  registerDoParallel(cores = 7)

  mod_net$caret_models <-
    map2(mod_globals$full_formulas,
         mod_net$tuning$opt_lambda,
         ~train(form = .x,
                data = stigma$analysis_tbl,
                method = 'glmnet',
                metric = 'RMSE',
                tuneGrid = data.frame(alpha = mod_net$alpha,
                                      lambda = .y),
                trControl = trainControl(method = 'repeatedcv',
                                         number = 10,
                                         savePredictions = 'final',
                                         returnData = TRUE,
                                         returnResamp = 'final'))) %>%
    map(as_caretx)

  stopImplicitCluster()

# Registering the futures ------

  insert_msg('Futures')

  plan('multisession')

# Model assumptions -------

  insert_msg('Model assumptions')

  ## plots of residuals

  mod_net$resid_plots <- mod_net$caret_models %>%
    future_map(plot,
               cust_theme = globals$common_theme,
               .options = furrr_options(seed = TRUE))

# Cross validation errors and R-squared ------

  insert_msg('Cross-validation fit stats')

  mod_net$cv_stats <- mod_net$caret_models %>%
    extract_caret_stats

# Plotting the R-squares and RMSE for the training and CV datasets ------

  insert_msg('Plotting the data and CV stats')

  mod_net$fit_stat_plots[c('RMSE', 'rsq')] <-
    plot_caret_stats(mod_net$cv_stats,
                     subtitle = 'Elastic Net') %>%
    map(~.x +
          scale_y_discrete(limits = rev(names(mod_net$caret_models)),
                           labels = exchange(names(mod_net$caret_models),
                                             dict = globals$stigma_lexicon,
                                             key = 'variable',
                                             value = 'label')))

# Non-zero model coefficients -----

  insert_msg('Non-zero model coefficients')

  ## extraction regex

  mod_net$extr_regex <- stigma$expl_lexicon$variable %>%
    sort(decreasing = TRUE) %>%
    paste(collapse = '|')

  ## coefficients

  mod_net$coeffs <- mod_net$glmnet_models %>%
    map(coef) %>%
    map(as.matrix) %>%
    map(~.x[.x != 0, ]) %>%
    map(as.data.frame,
        optional = TRUE,
        make.names = FALSE) %>%
    map(set_names, c('estimate')) %>%
    map(rownames_to_column, 'parameter') %>%
    map(mutate,
        variable = stri_extract(parameter,
                                regex = mod_net$extr_regex),
        level = stri_replace(parameter,
                             regex = mod_net$extr_regex,
                             replacement = ''),
        label = exchange(variable,
                         dict = stigma$expl_lexicon),
        ax_label = ifelse(parameter == '(Intercept)',
                          NA,
                          ifelse(level == '' | level == 'yes',
                                 label,
                                 paste(label, level, sep = ': ')))) %>%
    map(as_tibble)

  ## baselines

  mod_net$baselines <-  mod_net$coeffs %>%
    map(filter, parameter == '(Intercept)') %>%
    map_dbl(~.x$estimate[1])

# Plotting the coefficients ------

  insert_msg('Plotting the coefficients')

  mod_net$coef_plots <-
    list(x = mod_net$coeffs,
         y = exchange(names(mod_net$coeffs),
                      dict = globals$stigma_lexicon),
         z = mod_net$baselines) %>%
    pmap(function(x, y, z) x %>%
           filter(complete.cases(.)) %>%
           ggplot(aes(x = estimate,
                      y = reorder(ax_label, estimate),
                      size = abs(estimate),
                      color = factor(sign(estimate)))) +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_point(shape = 16) +
           geom_text(aes(label = signif(estimate, 3)),
                     size = 2.75,
                     vjust = -1.4,
                     show.legend = FALSE) +
           scale_color_manual(values = c('-1' = 'steelblue3',
                                         '1' = 'firebrick3'),
                              labels = c('-1' = 'negative',
                                         '1' = 'positive'),
                              name = expression(beta)) +
           scale_size_area(max_size = 4.5) +
           globals$common_theme +
           theme(axis.title.y = element_blank(),
                 axis.title.x = element_markdown()) +
           labs(title = y,
                subtitle = 'Elastic Net, non-zero estimates',
                x = paste('\u03B2<sub>ElasticNet</sub>, baseline =',
                           signif(z, 3))))

# END -----

  plan('sequential')

  insert_tail()
