# Launches the entire pipeline

  library(soucer)

  print(source_all(c('import.R',
                     'exploration.R',
                     'modeling.R',
                     'testing.R',
                     'incov_analyses.R',
                     'paper.R'),
                   message = TRUE,
                   crash = TRUE))
