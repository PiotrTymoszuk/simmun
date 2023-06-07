# This script contains project globals

# libraries ----

  library(plyr)
  library(tidyverse)
  library(stringi)

# Data container ------

  globals <- list()

# graphics -----

  globals$common_text <- element_text(size = 8,
                                      face = 'plain',
                                      color = 'black')

  globals$common_margin <- ggplot2::margin(t = 5,
                                           l = 4,
                                           r = 2,
                                           unit = 'mm')

  globals$common_theme <-
    theme_classic() +
    theme(axis.text = globals$common_text,
          axis.title = globals$common_text,
          plot.title = element_text(size = 8,
                                    face = 'bold'),
          plot.subtitle = globals$common_text,
          plot.tag = element_text(size = 8,
                                  face = 'plain',
                                  color = 'black',
                                  hjust = 0,
                                  vjust = 1),
          plot.tag.position = 'bottom',
          legend.text = globals$common_text,
          legend.title = globals$common_text,
          strip.text = globals$common_text,
          strip.background = element_rect(fill = 'gray95',
                                          color = 'gray80'),
          plot.margin = globals$common_margin,
          panel.grid.major = element_line(color = 'gray90'))

# Variable lexicon, the local cohort ------

  insert_msg('Variable lexicon, the local cohort')

  globals$stigma_lexicon <-
    list(age = c('age', 'years', NA),
         log_age = c('log age', 'years', NA),
         sqrt_age = c('sqrt age', 'years', NA),
         sex = c('sex', NA, NA),
         psych_comorb = c('mental disorder', NA, NA),
         somatic_comorb = c('physical disorder', NA, NA),
         infection = c('SARS-CoV2', NA, NA),
         bmi_class = c('body mass index', NA, NA),
         smoking = c('smoking', NA, NA),
         alcohol = c('alcohol', NA, NA),
         hads_anx_score = c('HADS, anxiety score', NA, NA),
         hads_dpr_score = c('HADS, depression score', NA, NA),
         hads_signs = c('anxiety/depression signs', NA, NA),
         pss_stress_score = c('PSS-4, mental stress score', NA, NA),
         log_pss_stress_score = c('log PSS-4, mental stress score', NA, NA),
         sqrt_pss_stress_score = c('sqrt PSS-4, mental stress score', NA, NA),
         phq_dpr_score = c('PHQ-9, depression score', NA, NA),
         log_phq_dpr_score = c('log PHQ-9, depression score', NA, NA),
         sqrt_phq_dpr_score = c('sqrt PHQ-9, depression score', NA, NA),
         study_group = c('study group', NA, NA),
         neo = c('NEO', 'nmol/L', 'inflammation'),
         log_neo = c('log NEO', 'nmol/L', 'inflammation'),
         sqrt_neo = c('sqrt NEO', 'nmol/L', 'inflammation'),
         il6 = c('IL6', 'ng/L', 'inflammation'),
         il6_class = c('IL6', NA, 'inflammation'),
         crp = c('CRP', 'mg/dL', 'inflammation'),
         crp_class = c('CRP', NA, 'inflammation'),
         nlr = c('NLR', NA, 'inflammation'),
         log_nlr = c('log NLR', NA, 'inflammation'),
         sqrt_nlr = c('sqrt NLR', NA, 'inflammation'),
         trp = c('TRP', 'µmol/L', 'metabolism'),
         log_trp = c('log TRP', 'µmol/L', 'metabolism'),
         sqrt_trp = c('sqrt TRP', 'µmol/L', 'metabolism'),
         kyn = c('KYN', 'µmol/L', 'metabolism'),
         log_kyn = c('log KYN', 'µmol/L', 'metabolism'),
         sqrt_kyn = c('sqrt KYN', 'µmol/L', 'metabolism'),
         phe = c('PHE', 'µmol/L', 'metabolism'),
         log_phe = c('log PHE', 'µmol/L', 'metabolism'),
         sqrt_phe = c('sqrt PHE', 'µmol/L', 'metabolism'),
         tyr = c('TYR', 'µmol/L', 'metabolism'),
         log_tyr = c('log TYR', 'µmol/L', 'metabolism'),
         sqrt_tyr = c('sqrt TYR', 'µmol/L', 'metabolism'),
         kyn_trp = c('KYN/TRP', NA, 'metabolism'),
         log_kyn_trp = c('log KYN/TRP', NA, 'metabolism'),
         sqrt_kyn_trp = c('sqrt KYN/TRP', NA, 'metabolism'),
         phe_tyr = c('PHE/TYR', NA, 'metabolism'),
         log_phe_tyr = c('log PHE/TYR', NA, 'metabolism'),
         sqrt_phe_tyr = c('sqrt PHE/TYR', NA, 'metabolism'),
         no = c('NO', 'µmol/L', 'inflammation'),
         log_no = c('log NO', 'µmol/L', 'inflammation'),
         sqrt_no = c('sqrt NO', 'µmol/L', 'inflammation'),
         neutro = c('neutro', '%', 'inflammation'),
         log_neutro = c('log neutro', '%', 'inflammation'),
         sqrt_neutro = c('sqrt neutro', '%', 'inflammation'),
         anti_rbd = c('anti-RBD IgG', 'AU', NA),
         log_anti_rbd = c('log anti-RBD IgG', 'AU', NA),
         anti_rbd_class = c('anti-RBD IgG', NA, NA),
         sqrt_anti_rbd = c('sqrt anti-RBD IgG', 'AU', NA)) %>%
    compress(names_to = 'variable',
             values_to = 'lab_list') %>%
    mutate(label = map_chr(lab_list, ~.x[[1]]),
           unit = map_chr(lab_list, ~.x[[2]]),
           var_group = map_chr(lab_list, ~.x[[3]]),
           axis_label = ifelse(!is.na(unit),
                               paste(label, unit, sep = ', '),
                               label))

  ## indicating the analysis responses

  globals$stigma_lexicon <- globals$stigma_lexicon %>%
    mutate(response = ifelse(stri_detect(variable,
                                         regex = '(age)|(sex)|(hads)|(pss)|(study_group)|(phq)'),
                             'no', 'yes'),
           base_label = stri_replace(label,
                                     regex = '(log\\s{1})|(sqrt\\s{1})',
                                     replacement = ''))

# Labels and names of the study groups, STIGMA -------

  insert_msg('Labels and colors for the study groups')

  globals$stigma_labs <- c('SARS-CoV-2.HADS+' = 'HADS+ CoV+',
                           'healthy.HADS+' = 'HADS+ CoV-',
                           'SARS-CoV-2.HADS-' = 'HADS- CoV+',
                           'healthy.HADS-' = 'HADS- CoV-')

  globals$stigma_colors <- c('SARS-CoV-2.HADS+' = 'coral4',
                             'healthy.HADS+' = 'gray60',
                             'SARS-CoV-2.HADS-' = 'coral2',
                             'healthy.HADS-' = 'darkolivegreen3')

# Names and labels of variables of interest in the INCOV cohort -----

  insert_msg('Variable names and labels, INCOV cohort')

  globals$incov_proteins <- c('IL6_INF',
                              'IL10_INF',
                              'TNF_INF',
                              'IFNG_INF')

  globals$incov_metabolites <- c('tryptophan',
                                 'kynurenine',
                                 'phenylalanine',
                                 'tyrosine',
                                 'dopamine 3-O-sulfate',
                                 'serotonin',
                                 'quinolinate')

# Labels for the time points, INCOV cohort ------

  insert_msg('Timepoint labels, INCOV cohort')

  globals$time_labs <- c('healthy' = 'uninfected',
                         'T1' = 'acute',
                         'T2' = 'sub-acute',
                         'T3' = 'recovery')

# END -----

  insert_tail()
