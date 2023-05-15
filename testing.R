# Univariable testing and correlation for factors found to be significantly
# associated with the TRP/KN and the PHE/TYR systems as indentified by
# linear modeling with backward elimination

# tools ------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(rlang)

  library(exda)

  library(soucer)

  insert_head()

  c('./tools/tools.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis globals ------

  insert_msg('Analysis globals')

  test_globals <- list()

  test_globals$variables <-
    list(trp = mod_eli$inference[c("trp", "kyn", "kyn_trp")],
         tyr = mod_eli$inference[c("tyr", "phe_tyr")]) %>%
    map(~map(.x,
             filter,
             p_value < 0.05,
             variable != 'Intercept')) %>%
    map(~map(.x, ~.x$variable)) %>%
    map(reduce, union)

# analysis scripts -----

  insert_msg('Analysis scripts')

  c('./testing scripts/trp_system.R',
    './testing scripts/tyr_system.R',
    './testing scripts/result_tables.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
