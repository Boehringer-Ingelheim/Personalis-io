# PersonalisIO

**This package provides convenience functions for reading real-world evidence data provided by [Personalis](https://www.personalis.com/)
into Bioconductor [MultiAssayExperiment](https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html) objects**

[![R-CMD-check](https://github.com/Boehringer-Ingelheim/Personalis-io/actions/workflows/test.yml/badge.svg)](https://github.com/Boehringer-Ingelheim/Personalis-io/actions/workflows/test.yml)
[![license](https://img.shields.io/badge/license-GPL3-blue.svg)](https://github.com/Boehringer-Ingelheim/Personalis-io/blob/master/LICENSE)
[![docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://Boehringer-Ingelheim.github.io/Personalis-io)

## Basic Usage

```r
# List of folders with Personalis data (each folder represents one sample)
sample_list = c("/tmp/personalis_example_data/sample1", "/tmp/personalis_example_data/sample2")

# Load the data
mae <- PersonalisIO::read_personalis(sample_list)

# mae is a MultiAssayExperiment with the following structure
mae
#> A MultiAssayExperiment object of 8 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 8:
#>  [1] gene_expression: SummarizedExperiment with 22956 rows and 5 columns
#>  [2] dna_small_variants_somatic: SummarizedExperiment with 6097 rows and 5 columns
#>  [3] dna_small_variants_tumor: SummarizedExperiment with 376654 rows and 5 columns
#>  [4] rna_small_variants_somatic: SummarizedExperiment with 5302 rows and 5 columns
#>  [5] rna_small_variants_tumor: SummarizedExperiment with 201347 rows and 5 columns
#>  [6] msi: SummarizedExperiment with 5 rows and 5 columns
#>  [7] cnv: SummarizedExperiment with 145 rows and 5 columns
#>  [8] tcr: SummarizedExperiment with 2 rows and 5 columns
```

For more details, see the [Documentation](https://boehringer-ingelheim.github.io/Personalis-io/):

-   [Vignette](https://boehringer-ingelheim.github.io/Personalis-io/articles/PersonalisIO.html)
-   [function reference](https://boehringer-ingelheim.github.io/Personalis-io/reference/index.html)

## Installation

For now, please install the development version of the package using `remotes`:

```R
install.packages("remotes")
remotes::install_github("Boehringer-Ingelheim/Personalis-io")
```

If there is enough interest, we will consider submitting the package to Bioconductor.

## Contact

Please use the [issue tracker](https://github.com/Boehringer-Ingelheim/Personalis-io/issues).

## Credits

This package was developed by [Gregor Sturm](https://github.com/grst) with contributions from [Christopher Mohr](https://github.com/christopher-mohr), Neetika Nath and Germán Leparc.
