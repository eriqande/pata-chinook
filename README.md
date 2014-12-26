# pata-chinook

A small repository with work that I am doing to test whether the somewhat
"off-label" use of GSI for identifying the origin of introduced populations
gives reasonable results when using our big chinook SNP baseline.  This is
in support of a project that Javier Ciancio is doing on the invasion of
Patagonian rivers by Chinook salmon.  

To reproduce this work.  Here is how you do it:

1. clone the repository
    ```sh
    git clone https://github.com/eriqande/pata-chinook.git
    ```
1. Open the RStudio project in the repo and first get the libraries that you need.
    You will need these:
    ```r
    library(digest)
    library(gpiper)
    library(stringr)
    library(reshape2)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    ```
    All of these are available from CRAN, _except for gpiper_.  To get that 
    package you must have the `devtools` package installed. Then issue this
    command:
    ```r
    devtools::install_github("eriqande/gpiper", ref = "3d5a6b5")
    ```
1. Once those libraries are installed, issue this at the R prompt:
    ```r
    source("R-main/gsi-off-label-check-sims.R")
    ```
1. The above will create the figure and put it where it needs to be.  To typeset
    the supplement, just run LaTeX/BibTeX on `tex/pata-chinook-suppl.tex`


__NOTE:__ The above three steps regenerate figures from cached results from the
simulations.  If you want to re-run the simulations then you have to modify the
line in `R-main/gsi-off-label-check-sims.R` that says:
```r
REDO_SIMS = FALSE
```
to:
```r
REDO_SIMS = TRUE
```
That will require a good several hours of computation...

## Terms 

As a work partially of the United States Government, this RStudio project is in the
public domain within the United States. Additionally, I waive
copyright and related rights in the work worldwide through the CC0 1.0
Universal public domain dedication.

The data included in this package are not works of the US Govt and are provided under
a different open source license that is yet to be determined.

See TERMS.md for more information.

