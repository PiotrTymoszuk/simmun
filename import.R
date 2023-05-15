# Import and wrangling of the manuscript data
#
# The local STIGMA cohort dataset encompasses readouts of inflammation
# and metabolic markers of neurotransmitter precursors along with the
# explanatory variables (anxiety/depression signs and persistent symptoms)
# and confounders (age and sex)
#
# The INCOV validation dataset contains selected metabolites of the TRP decay
# pathway (so called TRYCATS pathway) and key inflammatory cytokines measured
# in blood of healthy individuals and COVID-19 participants at consecutive time
# points (acute, sub-acute and recovery) after clnical onset.


# tools --------

  library(plyr)
  library(tidyverse)
  library(soucer)
  library(trafo)
  library(stringi)
  library(readxl)
  library(foreign)

  insert_head()

  c('./tools/globals.R',
    './tools/tools.R') %>%
    source_all(message = TRUE, crash = TRUE)

# Sourcing the import scripts for single studies ------

  insert_msg('Sourcing the import scripts')

  source_all(c('./import scripts/local.R',
               './import scripts/incov.R'),
             message = TRUE, crash = TRUE)

# Updating the project globals ------

  insert_msg('Project globals')

  globals$incov_lexicon <-
    rbind(incov$annotation_proteome %>%
            filter(variable %in% globals$incov_proteins),
          incov$annotation_metabolome %>%
            filter(variable %in% globals$incov_metabolites))

# SIMMUN, modeling and analysis dataset -----

  insert_msg('Modeling SIMUN dataset')

  ## modeling responses, the optimal transformations defined
  ## during explorative analysis

  stigma$responses <-
    c('trp' = 'trp',
      'kyn' = 'log_kyn',
      'kyn_trp' = 'log_kyn_trp',
      'phe' = 'log_phe',
      'tyr' = 'log_tyr',
      'phe_tyr' = 'sqrt_phe_tyr')

  ## explanatory variables

  stigma$expl_lexicon <-
    c(## infection-related features
      ## there are to few non-ambulatory COV cases
      ## for modeling
      'infection' = 'SARS-CoV-2',
      'anti_rbd_class' = 'anti-RBD IgG',
      #'severity' = 'COVID-19 severity',
      ## inflammatory markers
      ## IL6 and CRP are excludd from analysis
      ## they are elevated only in few participants
      ## total neutrophisl correlate nearly absolutely with NLR
      'log_neo' = 'log NEO, nmol/L',
      #'crp_class' = 'CRP',
      #'il6_class' = 'IL6',
      'log_nlr' = 'log NLR',
      #'sqrt_neutro' = 'log Neutro',
      ## demographics
      'age' = 'Age, decades',
      'sex' = 'Sex',
      ## psychometrics
      'pss_stress_score' = 'Stress, PSS-4',
      'hads_anx_score' = 'Anxiety, HADS',
      'hads_dpr_score' = 'Depression, HADS',
      ## comorbidity, medical history
      'psych_comorb' = 'Mental illness',
      'somatic_comorb' = 'Somatic illness',
      'bmi_class' = 'BMI',
      'smoking' = 'Smoking',
      'alcohol' = 'Alcohol') %>%
    compress(names_to = 'variable',
             values_to = 'label')

  ## analysis data table with complete data

  stigma$analysis_tbl <- stigma$data %>%
    mutate(infection = car::recode(cov,
                                   "'healthy' = 'no';
                                    'SARS-CoV-2' = 'yes'"),
           infection = factor(infection,
                              c('no', 'yes')),
           age = age/10) %>%
    select(patient_id,
           all_of(unname(stigma$responses)),
           all_of(stigma$expl_lexicon$variable)) %>%
    filter(complete.cases(.))

  ## complete IDs

  stigma$complete_ids <- stigma$analysis_tbl$patient_id

# INCOV: modeling and analysis dataset ------

  insert_msg('INCOV: modeling and analysis dataset')

  ## with complete cases

  incov$analysis_tbl <- incov[c("metabolome", "proteome")] %>%
    map(select,
        patient_id,
        timepoint,
        any_of(globals$incov_lexicon$variable)) %>%
    reduce(inner_join, by = c('patient_id', 'timepoint')) %>%
    right_join(incov$clinic[c('patient_id', 'timepoint', 'time_po',
                              'age', 'sex', 'bmi_class', 'bmi',
                              'cov', 'who_severity', 'severity')] %>%
                 mutate(time_po = ifelse(timepoint == 'healthy', 0, time_po)),
               .,
               by = c('patient_id', 'timepoint')) %>%
    filter(complete.cases(.))

  incov$complete_ids <- incov$analysis_tbl$patient_id

# END -----

  insert_tail()
