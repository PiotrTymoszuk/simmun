# Renders the manuscript and the supplementary material

  insert_head()

# bibliography ------

  insert_msg('Bibliography')

  cov_biblio <- read_bib('./paper/markdown/infl_biblio.bib')

# feedback and issues ------

  insert_msg('Rendering the issue file')

  render('./paper/markdown/issues.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# supplementary material ------

  insert_msg('Rendering the supplements')

  render('./paper/markdown/supplementary_material.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# paper -----

  insert_msg('Rendering the paper')

  render('./paper/markdown/main_manuscript.Rmd',
         output_format = my_word(),
         output_dir = './paper')

# rebuttal letter -----

  insert_msg('Rendering the Rebuttal Letter')

  render('./paper/markdown/rebuttal_letter.Rmd',
         output_format = word_document2(number_sections = FALSE,
                                        reference_docx = 'ms_template.docx'),
         output_dir = './paper')

# END -----

insert_tail()
