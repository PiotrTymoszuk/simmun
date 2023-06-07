# simmun
Metabolism of serum neurotransmitter precursors after COVID-19

<br>

## Summary

Effects of demographic, clinical, psychometric, inflammatory and infection-related variables on blood markers of systemic serotonin availability (serotonin, tryptophan), kynureine pathway products competing with serotonin biosynthesis (kynurenine, quinolinate) and on markers of systemic dopamine turnover (dopamine 3-O-sulfate, phenylalanine and tyrosine) were assessed in two cohorts of SARS-CoV-2 indidivuals and uninfected controls: the cross-sectional single-timepoint SIMMUN cohort and the published longitudinal INCOV collective [^1].

<br>

<p align = "center"> 
<img src = "https://github.com/PiotrTymoszuk/simmun/assets/80723424/b0ba4a3f-a36a-44fe-8ae9-13dc14112228" width = "80%">
</p>

Our results indicate that SARS-CoV-2-dependent inflammation can lower systemic availability of serotonin and dopamine by depletion of the tryptophan via the competitive kynurenine pathway and inhibition of the phenylalanine - tyrosine suppression, respectively. Those effects can be further amplified by advanced age, mental stress and depression. It remains to be investigated, if and how this mechanism may contribute to neurotransmitter metabolism in the central nervous system and, consequently, to psychiatric disorders following SARS-CoV-2 infection.

You may follow the analysis and manuscript progress [here](https://github.com/PiotrTymoszuk/simmun/tree/main/paper)

<br>

## Terms of use

Please cite the repository and the peer-reviewed publication, when available. The raw data files will be made upon request to the senior study author, [Prof. Katharina Hüfner](mailto:katharina.huefner@tirol-kliniken.at).

<br>

## Basic usage

The following development packages are required to run the pipeline:

```r

devtools::install_github('PiotrTymoszuk/soucer') ## script sourcing
devtools::install_github('PiotrTymoszuk/ExDA') ## exploratory data analysis and staristical hypothesis testing
devtools::install_github('PiotrTymoszuk/clustTools') ## factor analysis
devtools::install_github('PiotrTymoszuk/caretExtra') ## fit statistics and quality control for the Caret models
devtools::install_github('PiotrTymoszuk/lmqc') ## fit statistics and quality control for linear models
devtools::install_github('PiotrTymoszuk/figur') ## management of figures and tables in Rmd documents
devtools::install_github('PiotrTymoszuk/trafo') ## handling of tabular data

```

Source 'exec.R' to launch the entire pipeline:

```r

source('exec.R')

```
<br>

## Contact

The repository maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com). Data requests should be addressed to [Prof. Katharina Hüfner](mailto:katharina.huefner@tirol-kliniken.at).

<br>

## References

[^1]: Su, Y. et al. Multiple early factors anticipate post-acute COVID-19 sequelae. Cell 185, 881-895.e20 (2022).
