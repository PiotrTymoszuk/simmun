# Modeling of the TRP, KYN, KYN/TRP, PHE, TYR and PHE/TYR as a function
# of inflammatory, CoV-related and psychometric explanatory variables in
# the SIMMUN/STIGMA cohort
#
# Ramerkable over-fitting in LASSO!

# tools -------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(rlang)

  library(lmqc)
  library(MASS)
  library(glmnet)
  library(ranger)
  library(party)

  library(caret)
  library(caretExtra)

  library(exda)
  library(ggrepel)
  library(ggtext)

  library(doParallel)
  library(furrr)

  library(soucer)

  train <- caret::train
  explore <- exda::explore
  reduce <- purrr::reduce

  c('./tools/tools.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Analysis globals ------

  insert_msg('Analysis globals')

  mod_globals <- list()

  ## analysis table with normalized numeric variables

  mod_globals$analysis_tbl <- stigma$analysis_tbl %>%
    map_dfc(function(x) if(is.numeric(x)) scale(x)[, 1] else x)

  ## full and null model formulas

  mod_globals$full_formulas <- stigma$responses %>%
    map(~paste(.x,
               paste(stigma$expl_lexicon$variable, collapse = ' + '),
               sep = ' ~ ')) %>%
    map(as.formula)

  mod_globals$null_formulas <- stigma$responses %>%
    map(~paste(.x, '~ 1')) %>%
    map(as.formula)

# Analysis scripts ------

  insert_msg('Analysis scripts')

  c('./modeling scripts/elimination.R',
    './modeling scripts/glmnet.R',
    './modeling scripts/ranger.R',
    './modeling scripts/cforest.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END ------

  insert_tail()
