# Scripts to generate paper tables, figures and the supplement

# tools -----

  library(plyr)
  library(tidyverse)
  library(trafo)
  library(stringi)
  library(exda)

  library(cowplot)
  library(figur)
  library(rmarkdown)
  library(knitr)
  library(bookdown)
  library(flextable)

  library(soucer)

  source_all('./tools/tools.R',
             message = TRUE, crash = TRUE)

# paper scripts ------

  insert_msg('Sourcing the paper scripts')

  c('./paper scripts/tables.R',
    './paper scripts/figures.R',
    './paper scripts/supplement.R',
    './paper scripts/render.R') %>%
    source_all(message = TRUE, crash = TRUE)

# END -----

  insert_tail()
