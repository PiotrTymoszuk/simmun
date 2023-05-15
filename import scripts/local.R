

  insert_head()

# container -----

  stigma <- list()

# Reading the dataset from SPSS -----

  insert_msg('Reading the raw SPSS dataset')

  ## raw data for the selected variables

  stigma$raw <- read.spss('./data/local/01.06.2022.sav',
                          to.data.frame = TRUE)

  stigma$stress <- read.spss('./data/local/27.10.2021.sav',
                             to.data.frame = TRUE)

# Wrangling, NA handling, calculation of the ratios ------

  insert_msg('Wrangling')

  ## subscores for the stress and HADS depression/anxiety
  ## used later for testing of the psychometric tool consistency

  stigma$pss4 <- stigma$stress %>%
    mutate(patient_id = stri_replace(studien_id,
                                     regex = '\\s+$',
                                     replacement = '')) %>%
    select(patient_id,
           all_of(c('pss_4_1', 'pss_4_2', 'pss_4_3', 'pss_4_4')))

  stigma$pss4[c('pss_4_2', 'pss_4_3')] <-
    stigma$pss4[c('pss_4_2', 'pss_4_3')] %>%
    map_dfc(function(x) 4 - x)

  stigma[c('hads_depression',
           'hads_anxiety')] <-
    list(c('hads_depression_results',
           'hads_02', 'hads_04', 'hads_06', 'hads_08',
           'hads_10', 'hads_12', 'hads_14'),
         c('hads_anxiety_results',
           'hads_01', 'hads_03', 'hads_05', 'hads_07',
           'hads_09', 'hads_11', 'hads_13')) %>%
    map(~select(stigma$raw,
                studien_id,
                all_of(.x))) %>%
    map(mutate,
        patient_id = stri_replace(studien_id,
                                  regex = '\\s+$',
                                  replacement = '')) %>%
    map(select, -studien_id)

  ## stress, depression, final scores

  stigma$stress <- stigma$stress %>%
    as_tibble %>%
    transmute(patient_id = stri_replace(studien_id,
                                        regex = '\\s+$',
                                        replacement = ''),
              phq_dpr_score = phq_results,
              log_phq_dpr_score = log(phq_dpr_score + 1),
              sqrt_phq_dpr_score = sqrt(phq_dpr_score),
              pss_stress_score = pss_4_results,
              log_pss_stress_score = log(pss_stress_score + 1),
              sqrt_pss_stress_score = sqrt(pss_stress_score))


  ## inflammation , neurotransmitter precursors, demography, CoV
  ## antibodies

  stigma$data <- stigma$raw %>%
    as_tibble %>%
    transmute(patient_id = stri_replace(studien_id,
                                        regex = '\\s+$',
                                        replacement = ''),
              age = Alter_real,
              log_age = log(age),
              sqrt_age = sqrt(age),
              sex = car::recode(geschlecht,
                                "'Männlich' = 'male';
                                'Weiblich' = 'female'"),
              sex = factor(sex, c('female', 'male')),
              hads_anx_score = hads_anxiety_results,
              hads_dpr_score = hads_depression_results,
              hads_dpr_score = ifelse(is.na(hads_dpr_score),
                                      hads_anx_score,
                                      hads_dpr_score),
              hads_signs = ifelse(hads_anx_score >= 8 |
                                    hads_dpr_score >= 8,
                                  'HADS+', 'HADS-'),
              hads_signs = factor(hads_signs, c('HADS-', 'HADS+')),
              hads_anx_score = cut(hads_anx_score,
                                   c(-Inf, 7, Inf),
                                   c('< 8', '\u2265 8')),
              hads_dpr_score = cut(hads_dpr_score,
                                   c(-Inf, 7, Inf),
                                   c('< 8', '\u2265 8')),
              cov = car::recode(covid_real,
                                "'COVID pos' = 'SARS-CoV-2';
                                'COVID neg' = 'healthy'"),
              covid_normal_ward = ifelse(is.na(soziodemographische_daten_024_3),
                                         'no',
                                         as.character(soziodemographische_daten_024_3)),
              covid_normal_ward = car::recode(covid_normal_ward,
                                              "'ja' = 'yes'; 'nein' = 'no'"),
              covid_icu = ifelse(is.na(soziodemographische_daten_024_4),
                                 'no',
                                 as.character(soziodemographische_daten_024_4)),
              covid_icu = car::recode(covid_icu,
                                      "'ja' = 'yes'; 'nein' = 'no'"),
              severity = ifelse(cov == 'healthy',
                                'healthy',
                                ifelse(covid_icu == 'yes',
                                       'hospitalized',
                                       ifelse(covid_normal_ward == 'yes',
                                              'hospitalized', 'ambulatory'))),
              severity = factor(severity,
                                c('healthy', 'ambulatory', 'hospitalized')),
              psych_comorb = car::recode(psy_real,
                                         "'ja' = 'yes';
                                         'nein' = 'no'"),
              somatic_comorb = car::recode(aktuell_aktive_k_rperliche_erkrankungen,
                                           "'Ja' = 'yes'; 'nein' = 'no'"),
              study_group = interaction(cov, hads_signs),
              bmi_class = cut(bmi,
                              c(-Inf, 25, 30, Inf),
                              c('normal', 'overweight', 'obesity')),
              smoking = car::recode(nikotin, "'Ja' = 'yes'; 'Nein' = 'no'"),
              alcohol = car::recode(alkohol, "'Ja' = 'yes'; 'Nein' = 'no'"),
              new_psych_illness = ifelse(seit_welchem_jahr_besteht_psychiatr._hauptdiagnose__icd_. %in% c('2020', '2021'),
                                         'yes', 'no'),
              new_psych_illness = factor(new_psych_illness, c('no', 'yes')),
              neo = Neopterin,
              log_neo = log(neo),
              sqrt_neo = sqrt(neo),
              ## IL6 is stratified by the 7 ng/mL cutoff
              il6 = IL_6,
              il6_class = cut(il6,
                              c(-Inf, 7, Inf),
                              c('baseline', '> 7 pg/mL')),
              ## CRP is stratified by the 0.5 mg/dL cutoff
              crp = CRP,
              crp_class = cut(crp,
                              c(-Inf, 0.5, Inf),
                              c('baseline', '> 0.5 mg/dL')),
              nlr = NLR,
              log_nlr = log(nlr),
              sqrt_nlr = sqrt(nlr),
              trp = Tryptophan,
              log_trp = log(trp),
              sqrt_trp = sqrt(trp),
              kyn = Kynurenin,
              log_kyn = log(kyn),
              sqrt_kyn = sqrt(kyn),
              phe = Phenylalanin,
              log_phe = log(phe),
              sqrt_phe = sqrt(phe),
              tyr = Tyrosin,
              log_tyr = log(tyr),
              sqrt_tyr = sqrt(tyr),
              kyn_trp = kyn/trp,
              log_kyn_trp = log(kyn/trp),
              sqrt_kyn_trp = sqrt(kyn/trp),
              phe_tyr = phe/tyr,
              log_phe_tyr = log(phe/tyr),
              sqrt_phe_tyr = sqrt(phe/tyr),
              no = Nitrit,
              log_no = log(no),
              sqrt_no = sqrt(no),
              neutro = Segmentkernige_Neutrophile,
              log_neutro = log(neutro),
              sqrt_neutro = sqrt(neutro),
              anti_rbd = RBDpanIgIndex,
              log_anti_rbd = log(anti_rbd),
              sqrt_anti_rbd = sqrt(anti_rbd),
              ## stratification of the anti-RBD titer
              ## 0 - 1 negative (all healthy within the range)
              ## and by the median of titer in the CoV subset
              anti_rbd_class = cut(anti_rbd,
                                   c(-Inf, 1, 16.341369, Inf),
                                   c('negative', '1 - 16.3 AU', '> 16.3 AU')))

  stigma$data <- left_join(stigma$data,
                           stigma$stress,
                           by = 'patient_id')

# END ------

  stigma <- stigma[c('data', 'pss4', 'hads_depression', 'hads_anxiety')]

  rm(i)

  insert_tail()
