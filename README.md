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
# should take about 30 seconds
SKCM_stats = DCARSacrossNetwork(SKCM_rank,edgelist = STRING,
                                W = W, extractTestStatisticOnly = TRUE,
                                verbose = FALSE)
sort(SKCM_stats,decreasing=TRUE)[1:10]

# use speed-up to find significant (P = 0.05 unadjusted) gene pairs
# first take a stratified sample of gene pairs
sampleindices = stratifiedSample(SKCM_stats)
# then calculate their p-values using permutation testing
# this takes about 1 minute to run
# you can increase the number of iterations (niter) for a tighter graph
# but at the risk of longer runtimes
pvals = DCARSacrossNetwork(SKCM_rank,
                           edgelist = STRING[sampleindices,],
                           W = W, 
                           niter = 100,
                           verbose = FALSE)
# now calculate the critical value at P=0.05 using loess smoothing
criticalValue = getLoessCriticalValue(SKCM_stats[sampleindices], pvals, plot = TRUE)

# how many significant gene pairs do we have?
table(SKCM_stats > criticalValue)

# extract these significant edges
SKCM_signif_edges = STRING[SKCM_stats > criticalValue,]

# and graph them as a network
library(igraph)
SKCM_signif_graph = graph.edgelist(SKCM_signif_edges, directed = FALSE)
V(SKCM_signif_graph)$size = 0.01
V(SKCM_signif_graph)$label.cex = 0.4
plot(SKCM_signif_graph)

# extract network communities and plot
cm = walktrap.community(SKCM_signif_graph)
plot(cm, SKCM_signif_graph)
```

## Author

* **Shila Ghazanfar**  - [@shazanfar](https://twitter.com/shazanfar)

