# Handling of common links in the Rmarkdown project

  insert_head()

# container ------

  proj_links <- list()

# Links to development packages --------

  insert_msg('Links to development packages')

  proj_links <-
    list('ExDA' = 'https://github.com/PiotrTymoszuk/ExDA',
         'trafo' = 'https://github.com/PiotrTymoszuk/trafo',
         'figur' = 'https://github.com/PiotrTymoszuk/figur',
         'clustTools' = 'https://github.com/PiotrTymoszuk/clustTools',
         'lmqc' = 'https://github.com/PiotrTymoszuk/lmqc',
         'caretExtra' = 'https://github.com/PiotrTymoszuk/caretExtra') %>%
    compress(names_to = 'obj_name',
             values_to = 'x') %>%
    mutate(ref_name = paste0('_', obj_name, '_'))

  proj_links <- proj_links[c('ref_name', 'x')] %>%
    pmap(mdlink) %>%
    set_names(proj_links$obj_name)

# END ------

  insert_tail()
