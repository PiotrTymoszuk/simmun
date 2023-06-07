# Checking the consistency of the psychometric tools

  insert_head()

# container ----

  cons <- list()

# analysis globals -----

  insert_msg('Analysis globals')

  cons$analysis_tbl <-
    stigma[c("pss4", "hads_depression", "hads_anxiety")] %>%
    map(filter, patient_id  %in% stigma$complete_ids) %>%
    map(column_to_rownames, 'patient_id') %>%
    map(select, -ends_with('results')) %>%
    map(~filter(.x, complete.cases(.x)))

# Factor analysis -----

  insert_msg('Factor analysis')

  ## to check the tau equivalence

  cons$fa_objects <-
    list(data = cons$analysis_tbl,
         kdim = c(1, 3, 2)) %>%
    pmap(reduce_data,
         red_fun = 'fa')

  ## anxiety nearly tau equivalent! (all items in the vicinity of diagonal)

  cons$fa_loading_plots <-
    list(x = cons$fa_objects[c("hads_depression", "hads_anxiety")]) %>%
    pmap(plot,
         type = 'loadings',
         cust_theme = globals$common_theme) %>%
    map2(.,
         c('Depression, HADS',
           'Anxiety, HADS'),
         ~.x +
           labs(title = .y,
                subtitle = 'Factor analysis loadings') +
           geom_vline(xintercept = 0,
                      linetype = 'dashed') +
           geom_hline(yintercept = 0,
                      linetype = 'dashed') +
           geom_abline(slope = 1,
                       intercept = 0,
                       color = 'coral3',
                       linetype = 'dashed'))

# McDonald's omega ------

  insert_msg('McDonald omega')

  set.seed(1234)

  cons$omega_obj <-
    list(m = cons$analysis_tbl,
         nfactors = c(1, 3, 2)) %>%
    pmap(omega,
         fm = 'ml')

  cons$omega_stats <- cons$omega_obj %>%
    map(~.x$omega.group) %>%
    map(as.data.frame) %>%
    map(rownames_to_column, 'factor') %>%
    map(~.x[c('factor', 'total')]) %>%
    compress(names_to = 'tool') %>%
    as_tibble

# Plotting total omegas -----

  insert_msg('Plotting total omegas')

  cons$omega_plot <- cons$omega_stats %>%
    filter(factor == 'g') %>%
    ggplot(aes(x = total,
               y = reorder(tool, total))) +
    geom_bar(stat = 'identity',
             color = 'black',
             fill = 'steelblue') +
    geom_text(aes(label = signif(total, 2)),
              size = 2.75,
              hjust = 1.5,
              vjust = 0.5,
              color = 'white') +
    scale_y_discrete(labels = c(pss4 = 'PSS-4 mental stress',
                                hads_anxiety = 'HADS anxiety',
                                hads_depression = 'HADS depression')) +
    globals$common_theme +
    theme(axis.title.y = element_blank()) +
    labs(title = 'Psychometric tool consistency',
         subtitle = paste("Global McDonald's \u03C9, n =",
                          nrow(cons$analysis_tbl[[1]])),
         x = expression(omega))

# END ------

  insert_tail()
