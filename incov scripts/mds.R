# Visualizing cytokine and metabolite levels via multi-dimensional scaling
# Separate analyses for each timepoint and metabolite system
#
# In addition, correlations are visualized as a network plot

  insert_head()

# container ------

  incov_mds <- list()

# analysis tables ------

  insert_msg('Analysis tables')

  ## variables

  incov_mds$variables <-
    list(serotonin = c('serotonin',
                       'tryptophan',
                       'kynurenine',
                       'quinolinate',
                       'IL6_INF',
                       'IL10_INF',
                       'TNF_INF',
                       'IFNG_INF'),
         dopamine = c('dopamine.3.O.sulfate',
                      'phenylalanine',
                      'tyrosine',
                      'IL6_INF',
                      'IL10_INF',
                      'TNF_INF',
                      'IFNG_INF'))

  ## analysis tables with ranks

  incov_mds$analysis_tbl <- incov$analysis_tbl %>%
    select(patient_id,
           timepoint,
           all_of(globals$incov_lexicon$variable)) %>%
    set_names(make.names(names(.)))

  incov_mds$analysis_tbl <- incov_mds$variables %>%
    map(~select(incov_mds$analysis_tbl,
                patient_id,
                timepoint,
                all_of(.x))) %>%
    map(blast, timepoint) %>%
    unlist(recursive = FALSE) %>%
    map(select, -timepoint) %>%
    map(~map_dfc(.x, function(x) if(is.numeric(x)) rank(x) else x)) %>%
    map(column_to_rownames, 'patient_id') %>%
    map(t) %>%
    map(as.data.frame)

  ## variable lexicon

  incov_mds$lexicon <- globals$incov_lexicon %>%
    mutate(variable = make.names(variable))

# Calculation of distance - Spearman's correlation -----

  insert_msg('Scaled Spearman distance')

  ## is nothing else as a squared Euclidean distance between ranks
  ## scaling the results (min/max)

  incov_mds$dists <- incov_mds$analysis_tbl %>%
    map(calculate_dist, method = 'sumofsquares')

# MDS objects ------

  insert_msg('MDS objects')

  ## two dimensional MDS, wrapping into a red_analysis object

  incov_mds$mds_objects <-
    list(red_obj = incov_mds$dists %>%
           map(cmdscale, k = 2),
         component_tbl = incov_mds$dists %>%
           map(cmdscale, k = 2) %>%
           map(as.data.frame) %>%
           map(rownames_to_column, 'observation') %>%
           map(set_names, c('observation', 'comp_1', 'comp_2')),
         data = incov_mds$analysis_tbl %>%
           map(quo)) %>%
    transpose %>%
    map(~c(.x,
           list(loadings = NULL,
                red_fun = 'mds'))) %>%
    map(red_analysis)

# Plotting the MDS coordinates ------

  insert_msg('Plots of the coordinates')

  incov_mds$plots <-
    list(x = incov_mds$mds_objects,
         point_color = c(rep('indianred3', 4),
                         rep('darkolivegreen3', 4))) %>%
    pmap(plot,
         cust_theme = globals$common_theme)

  incov_mds$plots <-
    list(x = incov_mds$plots,
         y =  paste0(c('Uninfected', 'Acute Cov', 'Sub-acute CoV', 'Recovery'),
                     ', INCOV') %>%
           rep(2),
         z = paste('n =', map_dbl(incov_mds$analysis_tbl, ncol))) %>%
    pmap(function(x, y, z) x +
           labs(title = y,
                subtitle = z) +
           geom_text_repel(aes(label = exchange(observation,
                                                dict = incov_mds$lexicon)),
                           size = 2.5) +
           theme(plot.tag = element_blank()))

# Network analysis -------

  insert_msg('Network analysis')

  ## network objects created from correlation matrices
  ## (Spearman or Pearson on ranks in this case)
  ##
  ## weak effects (rho <= 0.2) are masked

  incov_mds$graph_obj <- incov_mds$analysis_tbl %>%
    map(t) %>%
    map(cor, method = 'pearson') %>%
    map(~ifelse(abs(.x) > 0.2, .x, 0)) %>%
    map(~(.x + 1)/2) %>%
    map(graph_from_adjacency_matrix,
        mode = 'undirected',
        diag = FALSE,
        weighted = TRUE)

# Force-directed network plots ------

  insert_msg('Force-directed network plots')

  ## colors

  incov_mds$net_colors <-
    ifelse(stri_detect(incov_mds$lexicon$variable,
                       regex = '_INF$'),
           'gray50',
           ifelse(incov_mds$lexicon$variable %in% incov_mds$variables$serotonin,
                  'orangered3', 'darkolivegreen4')) %>%
    set_names(incov_mds$lexicon$variable)

  ## plots

  set.seed(1234)

  incov_mds$network_plots <-
    list(graph = incov_mds$graph_obj,
         plot_title =  paste0(c('Uninfected', 'Acute Cov', 'Sub-acute CoV', 'Recovery'),
                              ', INCOV') %>%
           rep(2),
         plot_subtitle = paste('n =', map_dbl(incov_mds$analysis_tbl, ncol))) %>%
    pmap(draw_graph,
         show_edge_labs = FALSE)

# END -----

  insert_tail()
