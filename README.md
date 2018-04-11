DCARS <img src="man/figures/hex.png" align="right"  height="250" width="250"/>
======================================================

Differential Correlation across Ranked Samples

Overview
--------

**DCARS** is a flexible statistical approach which uses local weighted correlations to build a powerful and robust statistical test to identify significant variation in levels of concordance across a ranking of samples. This has the potential to discover biologically informative relationships between genes across a variable of interest, such as survival outcome.

Installation
--------

```r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("shazanfar/DCARS")
library(DCARS)
```

Usage
-----

The following example calculates the DCARS test statistic using the skin cutaneous melanoma (SKCM) TCGA dataset, along survival ranking. First it does so for a single gene pair, SKP1 and SKP2, and then EIF3C and EIF5B. Then DCARS selects the top 10 gene pairs after calculating statistics across the STRING network.

```r
data(STRING)
data(SKCM)
SKCM_rank = t(apply(SKCM,1,rank))

# highly significantly DCARS gene pair: SKP1 and SKP2
# calculates p-value based on permutation
DCARS(SKCM_rank,"SKP1","SKP2",plot=TRUE)
# extract only the test statistic
DCARS(SKCM_rank,"SKP1","SKP2", extractTestStatisticOnly = TRUE)

# not significantly DCARS gene pair: EIF3C and EIF5B
# calculates p-value based on permutation
DCARS(SKCM_rank,"EIF3C","EIF5B",plot=TRUE)
# extract only the test statistic
DCARS(SKCM_rank,"EIF3C","EIF5B", extractTestStatisticOnly = TRUE)

# build weight matrix
W = weightMatrix(ncol(SKCM_rank), type = "triangular", span = 0.5, plot = TRUE)

# extract DCARS test statistics
SKCM_stats = DCARSacrossNetwork(SKCM_rank,edgelist = STRING,
                                W = W, extractTestStatisticOnly = TRUE,
                                verbose = FALSE)
sort(SKCM_stats,decreasing=TRUE)[1:10]
```

## Author

* **Shila Ghazanfar**  - [@shazanfar](https://twitter.com/shazanfar)

