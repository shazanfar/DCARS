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
```

To calculate p-values associated with these tests quickly, we use permutation tests over a stratified sample of gene pairs

```r
# first take a stratified sample of gene pairs
sampleindices = stratifiedSample(SKCM_stats)

# What niter is needed for Bonferroni-corrected p-value of 0.05?
# Need to ensure that our choice of iterations is above this 'niter' value
# i do this by selecting 110% of niter
signif = 0.05
bonfSig = signif/length(SKCM_stats)
niter = (length(SKCM_stats)*(1/signif)) / length(sampleindices)
niter

# calculate a large number of permutation statistics from these stratified sample pairs
# this should take about 2 minutes
permstats = DCARSacrossNetwork(SKCM_rank,
                               edgelist = STRING[sampleindices,],
                               W = W, 
                               niter = round(1.1*niter),
                               verbose = TRUE,
                               extractPermutationTestStatistics = TRUE)
                               
permstatsVector = unlist(permstats)
permpvals = sapply(SKCM_stats, function(x) mean(permstatsVector >= x))
permpvals[permpvals==0] <- min(permpvals[permpvals!=0])/2

# plot DCARS statistic against estimated (unadjusted) p-value
# note if there are not enough iterations then the top-most may flatten out
plot(SKCM_stats, -log10(permpvals), col = "red", pch = 16)

# add bonferroni adjusted line
abline(h = -log10(bonfSig), lty = 2, col = "red")

# highlight FDR adjusted points
permpvals_FDR = p.adjust(permpvals, method = "BH")
FDR_sig = permpvals_FDR < 0.05
sum(permpvals_FDR < 0.05)

points(SKCM_stats[FDR_sig], -log10(permpvals)[FDR_sig], col = "blue", pch = 16, cex = 1.2)

# What are the most significant (FDR < 0.05) edges?
SKCM_signif_edges_FDR = STRING[FDR_sig,]
SKCM_signif_edges_FDR
```

Thresholding on unadjusted DCARS p-values (P-value < 0.05) to explore an unweighted network

```r
# extract these significant edges
SKCM_signif_edges = STRING[permpvals < 0.05,]

# and graph them as a network
library(igraph)
SKCM_signif_graph = graph.edgelist(SKCM_signif_edges, directed = FALSE)
V(SKCM_signif_graph)$size = 0.01
V(SKCM_signif_graph)$label.cex = 0.6
plot(SKCM_signif_graph)

# extract network communities and plot
cm = walktrap.community(SKCM_signif_graph)
plot(cm, SKCM_signif_graph)
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==1]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==2]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==3]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==4]))
```

Treating DCARS statistics as a network edge weight

```r
library(igraph)
SKCM_graph = graph.edgelist(STRING, directed = FALSE)
E(SKCM_graph)$weight = SKCM_stats

SKCM_meanStrength = strength(SKCM_graph, weights = E(SKCM_graph)$weight)/degree(SKCM_graph)

topNodes = names(degree(SKCM_graph))[order(SKCM_meanStrength, decreasing = TRUE)[1:20]]

plot(degree(SKCM_graph), SKCM_meanStrength)
text(degree(SKCM_graph)[topNodes], SKCM_meanStrength[topNodes],
     names(degree(SKCM_graph)[topNodes]))

library(ggplot2)
library(ggrepel)
df = data.frame(node = names(degree(SKCM_graph)),
                degree = degree(SKCM_graph), 
                meanStrength = SKCM_meanStrength, 
                topNodes = (names(degree(SKCM_graph)) %in% topNodes))
g = ggplot(df, aes(x = degree, y = meanStrength)) +
  geom_point() + 
  geom_text_repel(data = subset(df,topNodes==TRUE), aes(label = node)) +
  theme_minimal() +
  NULL
g

# this shows that the strongest genes are IMPDH1 and SKP2, these genes have the highest mean edge weights
```

## Author

* **Shila Ghazanfar**  - [@shazanfar](https://twitter.com/shazanfar)

