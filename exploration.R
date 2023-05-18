# Explorative data analysis
#
# characteristic of the cohorts
#
# Testing for the normality and EOV for the study variables (STIGMA and INCOV)
#
# Comparing included and excluded participants of the SIMMUN study

# tools --------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(stringi)
  library(exda)
  library(soucer)
  library(psych)
  library(clustTools)

  select <- dplyr::select
  reduce <- purrr::reduce

  insert_head()

  source_all('./tools/tools.R',
             message = TRUE, crash = TRUE)

# analysis scripts -------

  insert_msg('Analysis scripts')

  c('./exploration scripts/cohort.R',
    './exploration scripts/cohort_comparison.R',
    './exploration scripts/distribution.R',
    './exploration scripts/excluded.R',
    './exploration scripts/consistency.R',
    './exploration scripts/correlation.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
