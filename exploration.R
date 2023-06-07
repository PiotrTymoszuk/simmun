# Explorative data analysis
#
# characteristic of the cohorts
#
# Testing for the normality and EOV for the study variables (STIGMA and INCOV)
#
# Comparing included and excluded participants of the SIMMUN study
#
# Testing effects of SARS-CoV-2 and gender on metabolites
# and inflammatory markers (needed for the paper Discussion)

# tools --------

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(stringi)
  library(exda)
  library(soucer)
  library(psych)
  library(clustTools)
  library(lmqc)
  library(MASS)

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
    './exploration scripts/correlation.R',
    './exploration scripts/gender_simmun.R',
    './exploration scripts/gender_incov.R',
    './exploration scripts/severity_simmun.R',
    './exploration scripts/severity_incov.R',
    './exploration scripts/bmi_simmun.R',
    './exploration scripts/bmi_incov.R',
    './exploration scripts/age_simmun.R',
    './exploration scripts/age_incov.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
