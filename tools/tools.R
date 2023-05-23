# This script provides tools for project data analyses

# tools ----

  library(plyr)
  library(tidyverse)
  library(furrr)
  library(stringi)
  library(trafo)
  library(sciplot)
  library(rmarkdown)
  library(knitr)
  library(bookdown)
  library(clipr)

# Descriptive stats -----

  format_tbl <-function(data,
                        rm_complete = FALSE,
                        dict = NULL, ...) {

    data <- data %>%
      map_dfc(stri_replace,
              regex = '^no.*\\nyes:\\s{1}',
              replacement = '') %>%
      map_dfc(stri_replace,
              regex = '^Mean.*\\nMedian\\s{1}=\\s{1}',
              replacement = '') %>%
      map_dfc(stri_replace,
              fixed = 'Range',
              replacement = 'range') %>%
      map_dfc(stri_replace,
              fixed = 'Complete',
              replacement = 'complete')

    if(!is.null(dict)) {

      data <- data %>%
        exchange('variable',
                 dict = dict, ...)

    }

    if(rm_complete) {

      data <- data %>%
        map_dfc(stri_replace,
                regex = '\\ncomplete.*$',
                replacement = '')

    }

    data

  }

# modeling ------

  extract_caret_stats <- function(caret_list) {

    caret_list %>%
      future_map(summary.caretx,
                 .options = furrr_options(seed = TRUE)) %>%
      map(~map(.x, ~.x[c(3, 4), c('statistic', 'estimate')]) %>%
            compress(names_to = 'data_type')) %>%
      compress(names_to = 'response') %>%
      mutate(statistic = factor(statistic, c('RMSE', 'rsq')),
             data_type = factor(data_type, c('train', 'cv')))

  }

  plot_caret_stats <- function(caret_stats,
                               titles = c('Fit error, SIMMUN', 'Explained variance, SIMMUN'),
                               x_labs = list('RMSE', expression('R'^2)),
                               subtitle = NULL) {

    ## plots R-squared and RMSE for the training data and cross-validation

    list(data = blast(caret_stats, statistic),
         title = titles,
         x_lab = x_labs) %>%
      pmap(function(data, title, x_lab) data %>%
             ggplot(aes(x = estimate,
                        y = response,
                        fill = data_type)) +
             geom_bar(stat = 'identity',
                      position = position_dodge(0.9),
                      color = 'black') +
             # geom_text(aes(label = signif(estimate, 2)),
             #          size = 2.5,
             #         color = 'white',
             #        hjust = 1.3,
             #       position = position_dodge(0.9)) +
             scale_fill_manual(values = c(train = 'darkolivegreen4',
                                          cv = 'steelblue3'),
                               labels = c(train = 'training',
                                          cv = 'CV'),
                               name = '') +
             globals$common_theme +
             theme(axis.title.y = element_blank()) +
             labs(title = title,
                  x = x_lab,
                  subtitle = subtitle))

  }

# multiple testing -----

  re_adjust <- function(data, method = 'BH') {

    ## adjusts for multiple testing e.g. with the Benjamini-Hochberg method

    if(method != 'none') {

      data <- data %>%
        mutate(p_adjusted = p.adjust(p_value, method = method))

    }

    data %>%
      mutate(significance = ifelse(p_adjusted < 0.001,
                                   'p < 0.001',
                                   ifelse(p_adjusted >= 0.05,
                                          paste0('ns (p = ',
                                                 signif(p_adjusted, 2), ')'),
                                          paste('p =',
                                                signif(p_adjusted, 2)))))

  }

# Correlation bubble plot -------

  corr_buble <- function(data,
                         plot_title = NULL,
                         signif_only = TRUE,
                         rotate_x_labs = FALSE,
                         txt_size = 2.5,
                         fill_title = expression(rho)) {

    ## Bubble plot with Spearman's correlation coefficients

    data <- data %>%
      mutate(significant = ifelse(p_adjusted < 0.05, 'signif', 'ns'),
             fontface = ifelse(significant == 'signif', 'bold', 'plain'))

    if(signif_only) {

      data <- data %>%
        mutate(estimate = ifelse(p_adjusted < 0.05, estimate, NA))

      if(all(is.na(data[['estimate']]))) {

        warning('No significant correlations to plot.', call. = FALSE)

        return(NULL)

      }

    }

    plot <- data %>%
      ggplot(aes(x = variable1,
                 y = variable2,
                 size = abs(estimate),
                 fill = estimate)) +
      geom_point(shape = 21) +
      geom_text(aes(label = signif(estimate, 2),
                    alpha = significant,
                    fontface = fontface),
                size = txt_size,
                hjust = 0.5,
                vjust = -1.4) +
      scale_alpha_manual(values = c(ns = 0.45,
                                    signif = 1),
                         labels = c(ns = 'ns',
                                    signif = 'p < 0.05')) +
      scale_radius(limits = c(-1, 1),
                      range = c(0, 5)) +
      scale_fill_gradient2(low = 'steelblue',
                           mid = 'white',
                           high = 'firebrick',
                           limits = c(-1, 1),
                           name = fill_title) +
      scale_x_discrete(labels = set_names(globals$incov_lexicon$label,
                                          globals$incov_lexicon$variable)) +
      scale_y_discrete(labels = set_names(globals$incov_lexicon$label,
                                          globals$incov_lexicon$variable)) +
      guides(size = 'none',
             alpha = 'none') +
      globals$common_theme +
      theme(axis.title = element_blank()) +
      labs(title = plot_title,
           subtitle = paste('n =', data$n[1]))

    if(rotate_x_labs) {

      plot <- plot +
        theme(axis.text.x = element_text(hjust = 1, angle = 90))

    }

    plot

  }

# displaying post-hoc testing results: one-way plots ------

  add_p_one <- function(plot,
                        test,
                        txt_size = 2.5,
                        line_offset = 0.055,
                        txt_vjust = -0.4,
                        offset_factor = 2.2) {

    ## adds significant pairwise tests
    ## between healthy controls and times after CoV

    test <- test %>%
      mutate(plot_cap = ifelse(p.adj > 0.1,
                               'ns',
                               ifelse(p.adj >= 0.05,
                                      paste0('ns (p = ', signif(p.adj, 2), ')'),
                                      paste('p =', signif(p.adj, 2)))))

    y_diff <- diff(range(plot$data$variable))

    ## healthy - acute comparison

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'acute'))

    y_acute<- max(y_tbl$variable) +
      y_diff * line_offset

    if(test$plot_cap[1] != 'ns') {

      plot <- plot +
        annotate('segment',
                 x = 1,
                 xend = 2,
                 y = y_acute,
                 yend = y_acute,
                 size = 0.5) +
        annotate('text',
                 label = test$plot_cap[1],
                 x = 1.5,
                 y = y_acute,
                 vjust = txt_vjust,
                 size = txt_size)

    } else {

      y_acute <- -Inf

    }

    ## healthy - sub-acute comparison

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'sub-acute'))

    y_sub<- max(y_tbl$variable) +
      y_diff * line_offset

    if(test$plot_cap[2] != 'ns') {

      if(y_sub < y_acute + y_diff * line_offset * offset_factor) {

        y_sub <-  y_acute + y_diff* line_offset * offset_factor

      }

      plot <- plot +
        annotate('segment',
                 x = 1,
                 xend = 3,
                 y = y_sub,
                 yend = y_sub,
                 size = 0.5) +
        annotate('text',
                 label = test$plot_cap[2],
                 x = 2,
                 y = y_sub,
                 vjust = txt_vjust,
                 size = txt_size)

    } else {

      y_sub <- -Inf

    }

    ## healthy - recovery comparison

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'recovery'))

    y_reco <- max(y_tbl$variable) +
      y_diff * line_offset

    y_max <- max(y_acute, y_sub)

    if(test$plot_cap[3] != 'ns') {

      if(y_reco < y_max + y_diff * line_offset * offset_factor) {

        y_reco <-  y_max + y_diff * line_offset * offset_factor

      }

      plot <- plot +
        annotate('segment',
                 x = 1,
                 xend = 4,
                 y = y_reco,
                 yend = y_reco,
                 size = 0.5) +
        annotate('text',
                 label = test$plot_cap[3],
                 x = 2.5,
                 y = y_reco,
                 vjust = txt_vjust,
                 size = txt_size)

    }

    return(plot)

  }

  add_p_chain <- function(plot,
                          test,
                          txt_size = 2.5,
                          line_offset = 0.055,
                          txt_vjust = -0.4) {

    ## adds significant pairwise tests in a chain-wise manner
    ## group1 - group2, group2 - group3 etc.

    test <- test %>%
      mutate(plot_cap = ifelse(p > 0.1,
                               'ns',
                               ifelse(p.adj >= 0.05,
                                      paste0('ns (p = ', signif(p.adj, 2), ')'),
                                      paste('p =', signif(p.adj, 2)))))

    plot_var <- test[[1]][1]

    for(i in 1:nrow(test)) {

      if(test$plot_cap[i] == 'ns') {

        next

      }

      y_tbl <- plot$data %>%
        filter(group %in% c(test$group1[[i]],
                            test$group2[[i]]))

      y_line <- max(y_tbl$variable) +
        diff(range(y_tbl$variable)) * line_offset

      plot <- plot +
        annotate('segment',
                 x = i + 0.05,
                 xend = i + 1 - 0.05,
                 y = y_line,
                 yend = y_line,
                 size = 0.5) +
        annotate('text',
                 label = test$plot_cap[[i]],
                 x = i + 0.5,
                 y = y_line,
                 vjust = txt_vjust,
                 size = txt_size)


    }

    return(plot)

  }

  add_beta <- function(plot,
                       inference_data,
                       txt_size = 2.5,
                       line_offset = 0.055,
                       txt_vjust = -0.4,
                       offset_factor = 2.2,
                       p_sep = ', ') {

    ## adds model betas with confidence intervals and p values

    inference_data <- inference_data %>%
      mutate(est_lab = paste0(signif(estimate, 2), ' [',
                              signif(lower_ci, 2), ' - ',
                              signif(upper_ci, 2), ']'),
             plot_lab = paste(est_lab, significance, sep = p_sep))

    inference_labs <- inference_data %>%
      filter(variable != 'Intercept') %>%
      .$plot_lab

    ## appending the plots

    y_diff <- diff(range(plot$data$variable))

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'acute'))

    y_acute<- max(y_tbl$variable) +
      y_diff * line_offset

    plot <- plot +
      annotate('segment',
               x = 1,
               xend = 2,
               y = y_acute,
               yend = y_acute,
               size = 0.5) +
      annotate('text',
               label = inference_labs[1],
               x = 1.5,
               y = y_acute,
               vjust = txt_vjust,
               size = txt_size)

    ## healthy - sub-acute comparison

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'sub-acute'))

    y_sub <- max(y_tbl$variable) +
      y_diff * line_offset

    if(y_sub < y_acute + y_diff * line_offset * offset_factor) {

      y_sub <-  y_acute + y_diff* line_offset * offset_factor

    }

    plot <- plot +
      annotate('segment',
               x = 1,
               xend = 3,
               y = y_sub,
               yend = y_sub,
               size = 0.5) +
      annotate('text',
               label = inference_labs[2],
               x = 2,
               y = y_sub,
               vjust = txt_vjust,
               size = txt_size)

    ## healthy - recovery comparison

    y_tbl <- plot$data %>%
      filter(group %in% c('healthy', 'recovery'))

    y_reco <- max(y_tbl$variable) +
      y_diff * line_offset

    y_max <- max(y_acute, y_sub)

    if(y_reco < y_max + y_diff * line_offset * offset_factor) {

      y_reco <-  y_max + y_diff * line_offset * offset_factor

    }

    plot <- plot +
      annotate('segment',
               x = 1,
               xend = 4,
               y = y_reco,
               yend = y_reco,
               size = 0.5) +
      annotate('text',
               label = inference_labs[3],
               x = 2.5,
               y = y_reco,
               vjust = txt_vjust,
               size = txt_size)

    plot

  }

# Kinetic panels -----+

  plot_tc_ribbon <- function(data,
                             vars = c('IL6_INF', 'IL10_INF', 'TNF_INF', 'IFNG_INF'),
                             average_stat = c('mean', 'median'),
                             ribbon_stat = c('SEM', '2SEM', 'IQR'),
                             plot_title = NULL,
                             plot_subtitle = NULL,
                             y_lab = expression('normalized log'[2] * ' concentration'),
                             x_lab = 'SARS-CoV-2',
                             ribbon_alpha = 0.45) {

    ## plots average stats with errors as lines and ribbons

    ## entry ------

    average_stat <- match.arg(average_stat[1],
                              c('mean', 'median'))

    ribbon_stat <- match.arg(ribbon_stat[1],
                             c('SEM', '2SEM', 'IQR'))

    average_fun <- list(mean = function(x) mean(x, na.rm = TRUE),
                        median = function(x) median(x, na.rm = TRUE))

    ribbon_fun <- list(SEM = function(x) c(lower = mean(x, na.rm = TRUE) - se(x, na.rm = TRUE),
                                           upper = mean(x, na.rm = TRUE) + se(x, na.rm = TRUE)),
                       `2SEM` = function(x) c(lower = mean(x, na.rm = TRUE) - 2 * se(x, na.rm = TRUE),
                                              upper = mean(x, na.rm = TRUE) + 2 * se(x, na.rm = TRUE)),
                       IQR = function(x) c(lower = quantile(x, 0.25, na.rm = TRUE),
                                           upper = quantile(x, 0.75, na.rm = TRUE)))

    ## computing stats per timepoint -------

    data <- data[c('timepoint', vars)]

    summ_data <- data %>%
      blast(timepoint) %>%
      map(select, -timepoint) %>%
      map(~map(.x,
               ~c(stat = average_fun[[average_stat]](.x),
                  ribbon_fun[[ribbon_stat]](.x))) %>%
            reduce(rbind) %>%
            as.data.frame %>%
            mutate(variable = factor(vars, vars))) %>%
      compress(names_to = 'timepoint') %>%
      mutate(timepoint = factor(timepoint, levels(data[['timepoint']])))

    n_labs <- count(data, timepoint)

    n_labs <- map2_chr(n_labs[[1]],
                       n_labs[[2]],
                       paste, sep = '\nn = ') %>%
      set_names(n_labs[[1]])

    ## plotting -------

    summ_data %>%
      ggplot(aes(x = timepoint,
                 y = stat,
                 color = variable,
                 fill = variable)) +
      geom_hline(yintercept = 0,
                 linetype = 'dashed') +
      geom_path(aes(group = variable)) +
      geom_ribbon(aes(ymin = lower,
                      ymax = upper,
                      group = variable),
                  color = NA,
                  alpha = ribbon_alpha) +
      scale_x_discrete(labels = n_labs) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

  }

  plot_tc_beta <- function(data,
                           vars = c('IL6_INF', 'IL10_INF', 'TNF_INF', 'IFNG_INF'),
                           plot_title = NULL,
                           plot_subtitle = NULL,
                           y_lab = expression(beta * ', 95% CI'),
                           x_lab = 'SARS-CoV-2',
                           hide_baseline = TRUE,
                           baseline_lab = 'uninfected',
                           dodge_w = 0.5,
                           show_connector = FALSE,
                           connector_linetype = 'dashed',
                           show_estimates = FALSE,
                           txt_size = 2.5,
                           txt_hjust = 0.5,
                           txt_vjust = -0.6) {

    ## a Forest plot with modeling betas and 95% CI
    ## for the timepoints after the infection

    ## plot stats and data ------

    data <- data %>%
      map(mutate,
          level = ifelse(level == 'baseline',
                         baseline_lab, level)) %>%
      map(~mutate(.x,
                  level = factor(level, .x$level)))

    time_counts <- data[[1]][c('level', 'n', 'n_complete')]

    time_counts <- time_counts %>%
      mutate(n = ifelse(is.na(n),
                        time_counts$n_complete[1] - sum(n, na.rm = TRUE),
                        n),
             axis_lab = paste(level, n, sep = '\nn = '))

    n_labs <- set_names(time_counts$axis_lab,
                        time_counts$level)


    data <- data %>%
      reduce(rbind) %>%
      filter(response %in% vars) %>%
      mutate(response = factor(response, vars),
             est_lab = paste0(signif(estimate, 2), ' [',
                              signif(lower_ci, 2), ' - ',
                              signif(upper_ci, 2), ']'))

    if(hide_baseline) {

      data <- data %>%
        filter(level != baseline_lab)

    }

    ## plotting -------

    plot <- data %>%
      ggplot(aes(x = level,
                 y = estimate,
                 color = response)) +
      geom_hline(yintercept = 0,
                 linetype = 'dashed')

    if(show_connector) {

      plot <- plot +
        geom_path(aes(group = response),
                  linetype = connector_linetype,
                  position = position_dodge(width = dodge_w))

    }

    plot <- plot +
      geom_errorbar(aes(ymin = lower_ci,
                        ymax = upper_ci),
                    width = 0,
                    position = position_dodge(width = dodge_w)) +
      geom_point(shape = 16,
                 size = 2,
                 position = position_dodge(width = dodge_w)) +
      scale_x_discrete(labels = n_labs) +
      globals$common_theme +
      labs(title = plot_title,
           subtitle = plot_subtitle,
           x = x_lab,
           y = y_lab)

    if(show_estimates) {

      plot <- plot +
        geom_text(aes(label = est_lab),
                  size = txt_size,
                  hjust = txt_hjust,
                  vjust = txt_vjust,
                  angle = 90,
                  position = position_dodge(width = dodge_w))

    }

    plot

  }

# Visualization of adjacency ----

  draw_graph <- function(graph,
                         plot_title = NULL,
                         plot_subtitle = NULL,
                         show_edge_labs = TRUE,
                         rho_threshold = 0.2,
                         dict = incov_mds$lexicon,
                         node_colors = incov_mds$net_colors) {

    ## plots a force-directed graph

    net_plot <- graph %>%
      ggplot(aes(x = x,
                 y = y,
                 xend = xend,
                 yend = yend)) +
      geom_edges(aes(alpha = ifelse(abs(2*weight - 1) > rho_threshold,
                                    abs(2*weight - 1), 0),
                     size = ifelse(abs(2*weight - 1) > rho_threshold,
                                   abs(2*weight - 1), 0),
                     color = ifelse(abs(2*weight - 1) > rho_threshold,
                                    ifelse(weight > 0.5,
                                           'positive',
                                           ifelse(weight < 0.5,
                                                  'negative', 'ns')),
                                    'ns'))) +
      geom_nodes(aes(fill = name),
                 size = 3,
                 shape = 21,
                 color = 'white') +
      geom_nodelabel_repel(aes(label = exchange(name,
                                                dict = dict),
                               fill = name),
                           color = 'white',
                           size = 2.5,
                           show.legend = FALSE)

    if(show_edge_labs) {

      net_plot <- net_plot +
        geom_edgelabel_repel(aes(label = ifelse(abs(2*weight - 1) > rho_threshold,
                                                signif(2*weight - 1, 2),
                                                NA),
                                 color = ifelse(abs(2*weight - 1) > rho_threshold,
                                                ifelse(weight > 0.5,
                                                       'positive',
                                                       ifelse(weight < 0.5,
                                                              'negative', NA)),
                                                NA)),
                             size = 2.75,
                             show.legend = FALSE,
                             alpha = 1)

    }

    net_plot <- net_plot +
      scale_fill_manual(values = node_colors) +
      scale_size(limits = c(0, 1),
                 range = c(0.2, 1.7),
                 name = expression('abs(' * rho * ')')) +
      scale_alpha_continuous(limits = c(0, 1),
                             range = c(0.2, 1),
                             name = expression('abs(' * rho * ')')) +
      scale_color_manual(values = c(ns = 'gray90',
                                    positive = 'firebrick',
                                    negative = 'steelblue'),
                         name = 'Correlation') +
      guides(fill = 'none',
             color = guide_legend(),
             alpha = 'none') +
      theme_void() +
      theme(legend.text = globals$common_text,
            legend.title = globals$common_text,
            plot.subtitle = globals$common_text,
            plot.title = element_text(size = 8, face = 'bold')) +
      labs(title = plot_title,
           subtitle = plot_subtitle)


  }

# custom caret method for MM robust linear regression -----

  mm_rlm <-
    list(label = 'MM robust linear regression',
         library = 'MASS',
         type = 'Regression',
         parameters = data.frame(parameter = c('intercept', 'psi'),
                                 class = c('logical', 'character'),
                                 label = c('intercept', 'psi')),
         grid = function(x, y, search = "grid") {

           expand.grid(intercept = c(TRUE, FALSE), psi = c("psi.huber",
                                                           "psi.hampel", "psi.bisquare"))

         },
         fit = function (x, y, wts, param, lev, last, classProbs, ...) {
           dat <- if (is.data.frame(x))
             x
           else as.data.frame(x, stringsAsFactors = TRUE)
           dat$.outcome <- y
           psi <- MASS::psi.huber
           if (param$psi == "psi.bisquare")
             psi <- MASS::psi.bisquare
           else if (param$psi == "psi.hampel")
             psi <- MASS::psi.hampel
           if (!is.null(wts)) {
             if (param$intercept)
               out <- MASS::rlm(.outcome ~ ., data = dat, weights = wts,
                                psi = psi, method = 'MM', ...)
             else out <- MASS::rlm(.outcome ~ 0 + ., data = dat, weights = wts,
                                   psi = psi, method = 'MM', ...)
           }
           else {
             if (param$intercept)
               out <- MASS::rlm(.outcome ~ ., data = dat, psi = psi,
                                ...)
             else out <- MASS::rlm(.outcome ~ 0 + ., data = dat, psi = psi,
                                   ...)
           }
           out
         },
         predict = function (modelFit, newdata, submodels = NULL) {
           if (!is.data.frame(newdata))
             newdata <- as.data.frame(newdata, stringsAsFactors = TRUE)
           predict(modelFit, newdata)
         },
         prob = NULL,
         sor = function(x) x)

# markdown -----

  my_word <- function(...) {

    form <- word_document2(number_sections = FALSE,
                           reference_docx = 'ms_template.docx')

    form$pandoc$lua_filters <- c(form$pandoc$lua_filters,
                                 'scholarly-metadata.lua',
                                 'author-info-blocks.lua')

    form

  }

  insert_issue <- function(text = NULL) {

    if(!is.null(text)) {

      text <- paste0("<span custom-style = 'reviewer'>", text, "</span>")

    } else {

      text <- paste0("<span custom-style = 'reviewer'></span>")

    }

    write_clip(content = text,
               object_type = "character",
               breaks = "\n")

    return(text)

  }

  my_percent <- function(data,
                         variable,
                         return_n = FALSE,
                         signif_digits = 2) {

    ## percentages or count for categories of a qualitative variable

    counts <- table(data[[variable]])

    if(return_n) return(counts)

    complete <- sum(table(data[[variable]]))

    signif(counts/complete * 100,
           signif_digits)

  }

  my_beta <- function(data,
                      variables,
                      levels = NULL,
                      signif_digits = 2) {

    ## extracts beta with 95% confidence interval from an inference summary

    beta_dat <- data %>%
      filter(.data[['variable']] %in% variables)

    if(!is.null(levels)) {

      beta_dat <- beta_dat %>%
        filter(.data[['level']] == levels)

    }

    beta_dat %>%
      mutate(txt_lab = paste0(signif(estimate, signif_digits),
                              ' [95% CI: ', signif(lower_ci, signif_digits),
                              ' - ', signif(upper_ci, signif_digits), ']')) %>%
      .$txt_lab

  }

  my_rho <- function(data,
                     variables_x,
                     variables_y,
                     signif_digits = 2) {

    ## computes pair-wise Spearman's correlation coefficients

    cor_mtx <- data[, c(variables_x, variables_y)] %>%
      cor(method = 'spearman')

    cor_mtx[variables_x, variables_y] %>%
      signif(signif_digits)

  }

# END ------
