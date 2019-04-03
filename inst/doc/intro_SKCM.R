## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, warning  = FALSE,message = FALSE)
knitr::opts_chunk$set(cache = TRUE)

## ------------------------------------------------------------------------
data(STRING)
data(SKCM)
SKCM_rank = t(apply(SKCM,1,rank))

# highly significantly DCARS gene pair: SKP1 and SKP2
# calculates p-value based on permutation
DCARS(SKCM_rank,"SKP1","SKP2", plot=TRUE)

# extract only the test statistic
DCARS(SKCM_rank,"SKP1","SKP2", extractTestStatisticOnly = TRUE)
# examine the weighted correlation vector
wcor = DCARS(SKCM_rank,"SKP1","SKP2", extractWcorSequenceOnly = TRUE)
plot(wcor, type = "l", ylim = c(-1,1)); abline(h = 0, lty = 2)

# plot scatterplot split into five groups by survival ranking
plotColouredExpression(SKCM, genepair = c("SKP1","SKP2"), n = 5)

# plot ribbon plot of genes along sample ranking
plotOrderedExpression(SKCM, gene = c("SKP1", "SKP2"), facet = FALSE)


# not significantly DCARS gene pair: EIF3C and EIF5B
# calculates p-value based on permutation
DCARS(SKCM_rank,"EIF3C","EIF5B",plot=TRUE)
# extract only the test statistic
DCARS(SKCM_rank,"EIF3C","EIF5B", extractTestStatisticOnly = TRUE)
# examine the weighted correlation vector
wcor = DCARS(SKCM_rank,"EIF3C", "EIF5B",extractWcorSequenceOnly = TRUE)
plot(wcor, type = "l", ylim = c(-1,1)); abline(h = 0, lty = 2)

# build weight matrix
W = weightMatrix(ncol(SKCM_rank), type = "triangular", span = 0.5, plot = FALSE)

# extract DCARS test statistics
# should take about 30 seconds
SKCM_stats = DCARSacrossNetwork(SKCM_rank,edgelist = STRING,
                                W = W, extractTestStatisticOnly = TRUE,
                                verbose = FALSE)
sort(SKCM_stats,decreasing=TRUE)[1:10]

## ------------------------------------------------------------------------
globalCors = apply(STRING, 1, function(x) cor(SKCM_rank[x[1],], SKCM_rank[x[2],]))

# first take a stratified sample of gene pairs
sampleindices = stratifiedSample(abs(globalCors), length = 50)

# set number of permutations, roughly 1000 or so
niter = 1000

# calculate a large number of permutation statistics from these stratified sample pairs
# this should take about 3-4 minutes
permstats = DCARSacrossNetwork(SKCM_rank,
                               edgelist = STRING[sampleindices,],
                               W = W, 
                               niter = niter,
                               weightedConcordanceFunction = weightedPearson_matrix,
                               weightedConcordanceFunctionW = "matrix",
                               verbose = FALSE,
                               extractPermutationTestStatistics = TRUE)

res = estimatePvaluesSpearman(stats = SKCM_stats, 
globalCors = globalCors, permstats = permstats, 
usenperm = TRUE, nperm = 5000, plot = TRUE, verbose = FALSE)
                               
# plot DCARS statistic against estimated (unadjusted) p-value
# note if there are not enough iterations then the top-most may flatten out
plot(SKCM_stats, -log10(res$pval), col = "red", pch = 16)

# highlight FDR adjusted points
pval_FDR = p.adjust(res$pval, method = "BH")
FDR_sig = pval_FDR < 0.3
sum(pval_FDR < 0.3)

points(SKCM_stats[FDR_sig], -log10(res$pval)[FDR_sig], col = "blue", pch = 16, cex = 1.2)

# What are the most significant (FDR < 0.3) edges?
SKCM_signif_edges_FDR = STRING[FDR_sig,]
SKCM_signif_edges_FDR

## ------------------------------------------------------------------------
# extract these significant edges
SKCM_signif_edges = STRING[res$pval < 0.05,]

# and graph them as a network
library(igraph)
SKCM_signif_graph = graph.edgelist(SKCM_signif_edges, directed = FALSE)
V(SKCM_signif_graph)$size = 0.01
V(SKCM_signif_graph)$label.cex = 0.6
plot(SKCM_signif_graph)

# plot network with pathway information overlaid
# pathway information downloaded from MSigDB REACTOME pathways
plotNetworkPathway(SKCM_signif_graph)

# extract network communities and plot
cm = walktrap.community(SKCM_signif_graph)
plot(cm, SKCM_signif_graph)
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==1]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==2]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==3]))
plot(induced_subgraph(SKCM_signif_graph, V(SKCM_signif_graph)[cm$membership==4]))

## ------------------------------------------------------------------------
library(igraph)
SKCM_graph = graph.edgelist(STRING, directed = FALSE)
E(SKCM_graph)$weight = SKCM_stats
E(SKCM_graph)$pvalweight = -log10(res$pval)

# Graph the Ego network associated with SKP1 with 
plotEgoNetwork(hubnode = c("SKP1"), network = SKCM_graph, weight = "weight")

# only display unadjusted P < 0.05 significant pairs
plotEgoNetwork(hubnode = c("SKP1"), network = SKCM_graph, weight = "pvalweight", subset = TRUE, thresh = -log10(0.05))

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

## ------------------------------------------------------------------------
# calculate the weighted correlation vectors across the network
SKCM_wcor = t(DCARSacrossNetwork(SKCM_rank,
                               edgelist = STRING,
                               W = W, 
                               verbose = FALSE,
                               extractWcorSequenceOnly = TRUE))
plotWCorLine(SKCM_wcor, gene = "SKP1")

## ------------------------------------------------------------------------
SKCM_signif_wcor = t(DCARSacrossNetwork(SKCM_rank,
                                        edgelist = SKCM_signif_edges,
                                        W = W, 
                                        verbose = FALSE,
                                        extractWcorSequenceOnly = TRUE))

plotWcorsClusterPathway(SKCM_signif_wcor,cluster = TRUE, cutk = 6)

## ------------------------------------------------------------------------
sessionInfo()

