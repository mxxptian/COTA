# COTA

`COTA` is an R package to identify candidate core disease genes ssing trans-regulatory effects.

# Installation


You can install the development version of
`COTA` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

Before installation of COTA, you are also requested the below packages:
``` r
install.packages(c('qvalue', 'data.table', 'stringr', 'tidyr', 'AnnotationDbi', 'org.Hs.eg.db', 'ggplot2', 'igraph', 'VennDiagram', 'biomaRt', 'plyr', 'dplyr'), dependencies=TRUE)

```

``` r
devtools::install_github("mxxptian/COTA")
```


