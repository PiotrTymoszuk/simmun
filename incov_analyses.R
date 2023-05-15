# analyses in the INCOV cohort

# tools ------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(rlang)

  library(soucer)
  library(exda)
  library(ggpubr)

  library(MASS)
  library(lmqc)
  library(caret)
  library(clustTools)

  library(doParallel)
  library(furrr)

  library(igraph)
  library(ggnetwork)

  library(ggrepel)

  library(soucer)

  select <- dplyr::select
  reduce <- purrr::reduce
  explore <- exda::explore

  insert_head()

  c('./tools/tools.R') %>%
    source_all(message = TRUE, crash = TRUE)

# analysis scripts ------

  insert_msg('Analysis scripts')

  ## modeling of serum levels of dopamine sulfate and serotonin

  if(file.exists('./cache/incov_neuro.RData')) {

    insert_msg('Loading cached INCOV modeling results for serotonin and dopamine')

    load('./cache/incov_neuro.RData')

  } else {

    source_all('./incov scripts/neuro_modeling.R')

  }

  ## modeling for particular metabolites

  if(file.exists('./cache/incov_mod.RData')) {

    insert_msg('Loading cached INCOV modeling results')

    load('./cache/incov_mod.RData')

  } else {

    source_all('./incov scripts/modeling.R',
               message = TRUE, crash = TRUE)

  }

  ## univariable analyses

  c('./incov scripts/metabolite_correlation.R',
    './incov scripts/cytokine_correlation.R',
    './incov scripts/correlation.R',
    './incov scripts/mds.R',
    './incov scripts/time_course.R',
    './incov scripts/time_modeling.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
