# Import of the INCOV cohort proteomic and metabolomic data

  insert_head()

# container ----

  incov <- list()

# importing the metabolomics and proteomics data ------

  insert_msg('Metabolomics and proteomics data')

  incov[c('proteome',
          'metabolome')] <- c('S2.1 Proteomics',
                              'S2.2 Metabolomics') %>%
    map(~read_excel(path = './data/INCOV/Table S2.xlsx',
                    sheet = .x))

  incov[c('proteome',
          'metabolome')] <- incov[c('proteome',
                                    'metabolome')] %>%
    map(mutate,
        patient_id = `Patient Subject ID`,
        sample_id = ifelse(`Healthy or INCOV` == 'Healthy',
                           patient_id,
                           `Blood Draw`),
        .before = `Patient Subject ID`) %>%
    map(select,
        - `Patient Subject ID`,
        - `Blood Draw`,
        - `Healthy or INCOV`,
        - age,
        - sex,
        - BMI)

# Annotation tables -----

  insert_msg('Annotation tables')

  incov$annotation_proteome <-
    tibble(variable = names(incov$proteome)[!names(incov$proteome) %in% c('patient_id', 'sample_id')]) %>%
    mutate(label = stri_replace(variable, regex = '_.*$', replacement = ''),
           unit = 'Z-score of log\u2082 concentration',
           axis_lab = paste(label, unit, sep = ', '))

  incov$annotation_metabolome <-
    tibble(variable = names(incov$metabolome)[!names(incov$metabolome) %in% c('patient_id', 'sample_id')]) %>%
    mutate(label = variable,
           label = car::recode(label,
                               "'tryptophan' = 'TRP';
                               'kynurenine' = 'KYN';
                               'phenylalanine' = 'PHE';
                               'tyrosine' = 'TYR';
                               'dopamine 3-O-sulfate' = 'DA sulfate';
                               'serotonin' = 'serotonin';
                               'quinolinate' = 'QUIN'"),
           unit = 'Z-score of log\u2082 concentration',
           axis_lab = paste(label, unit, sep = ', '))

# Reading the clinical data for CoV patients -------

  insert_msg('Reading the clinical data for CoV patients')

  ## warnings concern coercion during reading from Excel
  ## and do not affect the data

  incov$clinic$cov <-
    suppressWarnings(read_excel('./data/INCOV/Table S1.xlsx',
                                sheet = 'S1.1 INCOV')) %>%
    transmute(patient_id = `Study subject ID`,
              sample_id = `Blood draw`,
              timepoint = stri_extract(sample_id, regex = 'T\\d{1}$'),
              time_po = `Observation days since onset of symptoms`,
              who_severity = `Who Ordinal Scale`,
              bmi = `BMI at baseline`,
              sex = tolower(Sex),
              sex = factor(sex, c('female', 'male')),
              age = `Age at baseline`,
              race = Race,
              cov = 'SARS-CoV-2')

  ## acute COVID-19 severity

  incov$clinic$acute <- incov$clinic$cov %>%
    filter(timepoint == 'T1') %>%
    mutate(acute_who_severity = who_severity,
           severity = car::recode(as.character(acute_who_severity),
                                  "'0' = 'mild';
                                  '1' = 'mild';
                                  '1 or 2' = 'mild';
                                  '3' = 'moderate';
                                  '4' = 'moderate';
                                  '5' = 'severe';
                                  '6' = 'severe';
                                  '7' = 'critical'"),
           severity = ifelse(severity == '<=2',
                             'mild', severity),
           bi_severity = ifelse(severity %in% c('mild', 'moderate'),
                                'mild/moderate',
                                ifelse(severity %in% c('severe', 'critical'),
                                       'severe/critical', NA))) %>%
    select(patient_id,
           acute_who_severity,
           severity,
           bi_severity)

  incov$clinic$cov <- left_join(incov$clinic$cov,
                                incov$clinic$acute,
                                by = 'patient_id')

  ## persistent symptom data

  incov$clinic$pasc <- read_excel('./data/INCOV/Table S1.xlsx',
                                  sheet = 'S1.3 PASC data') %>%
    mutate(patient_id = `Study subject ID`) %>%
    select(- `Study subject ID`)

  incov$clinic$pasc_present <- incov$clinic$pasc %>%
    select(- patient_id) %>%
    map_dfc(~ifelse(is.na(.x), 0, .x)) %>%
    reduce(`+`)

  incov$clinic$pasc <- incov$clinic$pasc %>%
    map_dfc(car::recode, "'1' = 'yes'; '0' = 'no'") %>%
    mutate(long_cov = ifelse(incov$clinic$pasc_present == 0, 'no', 'yes'),
           psy_long_cov = ifelse(depression == 'yes' |
                                   anxiety == 'yes' |
                                   difficulty_sleeping == 'yes',
                                 'yes', 'no'),
           psy_long_cov = ifelse(is.na(psy_long_cov),
                                 'no', psy_long_cov),
           neuro_long_cov = ifelse(Neurological == 'yes' |
                                     `Anosmia/Dysgeusia` == 'yes' |
                                     loss_of_taste == 'yes' |
                                     loss_of_smell == 'yes' |
                                     blurry_vision == 'yes' |
                                     pain_feet_hands == 'yes' |
                                     dizziness == 'yes',
                                   'yes', 'no'),
           neuro_long_cov = ifelse(is.na(neuro_long_cov),
                                   'no', neuro_long_cov),
           gastro_long_cov = ifelse(Gastrointestinal == 'yes' |
                                      diarrhea == 'yes' |
                                      nausea == 'yes',
                                    'yes', 'no'),
           gastro_long_cov = ifelse(is.na(gastro_long_cov),
                                    'no', gastro_long_cov),
           cogn_long_cov = ifelse(memory_problems == 'yes' |
                                    difficulty_concentrating == 'yes',
                                  'yes', 'no'),
           cogn_long_cov = ifelse(is.na(cogn_long_cov), 'no', cogn_long_cov),
           resp_long_cov = ifelse(`Respiratory viral` == 'yes' |
                                    cough == 'yes' |
                                    shortness_of_breath == 'yes',
                                  'yes', 'no'),
           resp_long_cov = ifelse(is.na(resp_long_cov),
                                  'no', resp_long_cov),
           fatigue_long_cov = ifelse(fatigue == 'yes' |
                                       inability_to_exercise == 'yes',
                                     'yes', 'no'),
           fatigue_long_cov = ifelse(is.na(fatigue_long_cov),
                                     'no', fatigue_long_cov))

  incov$clinic$cov <- left_join(incov$clinic$cov,
                                incov$clinic$pasc,
                                by = 'patient_id') %>%
  mutate(psy_long_cov = ifelse(psy_long_cov == 'no' & long_cov == 'no',
                               'recovered', psy_long_cov),
         neuro_long_cov = ifelse(neuro_long_cov == 'no' & long_cov == 'no',
                                 'recovered', neuro_long_cov),
         gastro_long_cov = ifelse(gastro_long_cov == 'no' & long_cov == 'no',
                                  'recovered', gastro_long_cov),
         cogn_long_cov = ifelse(cogn_long_cov == 'no' & long_cov == 'no',
                                'recovered', cogn_long_cov),
         resp_long_cov = ifelse(resp_long_cov == 'no' & long_cov == 'no',
                                'recovered', resp_long_cov),
         fatigue_long_cov = ifelse(fatigue_long_cov == 'no' & long_cov == 'no',
                                   'recovered', fatigue_long_cov))

  incov$clinic$pasc <- NULL
  incov$clinic$acute <- NULL
  incov$clinic$pasc_present <- NULL

  incov <- compact(incov)

# reading the clinical information for healthy individuals ------

  insert_msg('Healthy participant information')

  incov$clinic$healthy <-
    read_excel('./data/INCOV/Table S1.xlsx',
               sheet = 'S1.2 Healthy control') %>%
    transmute(patient_id = `Study Subject ID`,
              sample_id = patient_id,
              timepoint = 'healthy',
              bmi = BMI,
              sex = tolower(Sex),
              sex = factor(sex, c('female', 'male')),
              age = Age,
              race = Race,
              cov = 'healthy',
              who_severity = '0',
              acute_who_severity = '0',
              severity = 'healthy',
              long_cov = 'healthy',
              psy_long_cov = 'healthy',
              gastro_long_cov = 'healthy',
              neuro_long_cov = 'healthy',
              cogn_long_cov = 'healthy',
              resp_long_cov = 'healthy',
              fatigue_long_cov = 'healthy')

# Merging the healthy and CoV clinical information ------

  insert_msg('Merging the clinical information')

  incov$clinic <- full_rbind(incov$clinic$cov,
                             incov$clinic$healthy) %>%
    mutate(timepoint = car::recode(timepoint,
                                   "'healthy' = 'healthy';
                                   'T1' = 'acute';
                                   'T2' = 'sub-acute';
                                   'T3' = 'recovery'"),
           timepoint = factor(timepoint,
                             c('healthy', 'acute', 'sub-acute', 'recovery')),
           cov = factor(cov, c('healthy', 'SARS-CoV-2')),
           who_severity = factor(who_severity,
                                 c('0', '1', '2', '1 or 2',
                                   '<=2', '3', '4', '5', '6', '7')),
           acute_who_severity = factor(acute_who_severity,
                                       c('0', '1', '2', '1 or 2',
                                         '<=2', '3', '4', '5', '6', '7')),
           severity = factor(severity,
                             c('healthy', 'mild',
                               'moderate', 'severe', 'critical')),
           bi_severity = factor(bi_severity,
                                c('healthy',
                                  'mild/moderate',
                                  'severe/critical')),
           long_cov = factor(long_cov,
                             c('healthy', 'no', 'yes')),
           psy_long_cov = factor(psy_long_cov,
                                 c('healthy', 'recovered', 'no', 'yes')),
           neuro_long_cov = factor(neuro_long_cov,
                                   c('healthy', 'recovered', 'no', 'yes')),
           gastro_long_cov = factor(gastro_long_cov,
                                    c('healthy', 'recovered', 'no', 'yes')),
           cogn_long_cov = factor(cogn_long_cov,
                                  c('healthy', 'recovered', 'no', 'yes')),
           resp_long_cov = factor(resp_long_cov,
                                  c('healthy', 'recovered', 'no', 'yes')),
           fatigue_long_cov = factor(fatigue_long_cov,
                                     c('healthy', 'recovered', 'no', 'yes')),
           age = as.numeric(age),
           bmi = as.numeric(bmi),
           bmi_class = cut(bmi,
                           c(-Inf, 25, 30, Inf),
                           c('normal', 'overweight', 'obesity')),
           ethnics = car::recode(race,
                                 "'american indian or alaska native' = 'Other';
                                 'American Indian/Alaska Native' = 'Other';
                                 'ashkenazi jewish' = 'other';
                                 'asian' = 'Asian';
                                 'black or african-american' = 'Black or African-American';
                                 'Black or African American' = 'Black or African-American';
                                 'east asian' = 'Asian';
                                 'hispanic latino or spanish origin' = 'Hispanic';
                                 'Korean' = 'Asian';
                                 'middle eastern or north african' = 'Other';
                                 'More Than One Race' = 'Other';
                                 NA = 'Other';
                                 'Native Hawaiian or Other Pacific Islander' = 'Other';
                                 'Non Hispanic/Latino' = 'Other';
                                 'other' = 'Other';
                                 'south asian' = 'Asian';
                                 'Unknown / Not Reported' = 'Other';
                                 'Unreported' = 'Other';
                                 'white' = 'White';
                                 'White' = 'White'"),
           ethnics = factor(ethnics,
                            c('Asian',
                              'Black or African-American',
                              'White',
                              'Other')))

# appanding the omisc data with the clinical information ----

  insert_msg('Appending the omics data with clinical information')

  incov[c('proteome', 'metabolome')] <- incov[c('proteome', 'metabolome')] %>%
    map(~inner_join(incov$clinic, .x,
                    by = c('patient_id', 'sample_id')))

# END ------

  insert_tail()
