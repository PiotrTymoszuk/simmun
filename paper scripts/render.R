# Renders the manuscript and the supplementary material

  insert_head()

# bibliography ------

  insert_msg('Bibliography')

  cov_biblio <- read_bib('./paper/markdown/infl_biblio.bib')

# supplementary material ------

  insert_msg('Rendering the supplements')

  render('./paper/markdown/supplementary_material.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# paper -----

  insert_msg('Rendering the paper')

  render('./paper/markdown/main_manuscript.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# END -----

insert_tail()