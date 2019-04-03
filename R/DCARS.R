#' the weightMatrix function
#'
#' @title weightMatrix
#' @param n the number of samples, same as number of columns of data
#' @param type one of "harmonic", "triangular" or "block"
#' @param span proportion of samples to include on either side, default is 0.5
#' @param plot if TRUE, a heatmap of weight matrix will be generated
#' @return \code{matrix} an n*n matrix is generated corresponding to the weights for each sample

#' @examples
#'
#' weightMatrix(100)
#' weightMatrix(100, plot = TRUE)
#' weightMatrix(100, type = "triangular", span = 0.1, plot = TRUE)
#' weightMatrix(100, type = "block", plot = TRUE)
#' weightMatrix(100, type = "harmonic", plot = TRUE)
#'
#' @export

weightMatrix = function(n,
                        type = "triangular",
                        span = NULL,
                        plot = FALSE) {

  # for the number of samples, a n*n matrix is generated
  # corresponding to the weights for each sample
  # n is the number of samples
  # type can be one of "harmonic", "triangular" or "block"
  # default type is "triangular"
  # span is the proportion of samples to include on either side
  # default is 0.5

  if (is.null(type)) {
    type = "triangular"
    message("type not specified, defaulting to triangular")
  }

  if (!type%in%(c("harmonic","triangular","block"))) {
    type = "triangular"
    message("type given is not one of \"harmonic\",\"triangular\", or \"block\", defaulting to \"triangular\"")
  }

  if (type!="harmonic") {
  if (is.null(span)) {
    span = 0.5
    message("span not specified, defaulting to 0.5")
  }

  if (span<0 | span>1) {
    span = 0.5
    message("span specified not between 0 and 1, defaulting to 0.5")
  }
  }

  W = matrix(0,nrow=n,ncol=n)
  spanSamples = ceiling(span*n)

  if (type=="harmonic") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i,j] <- 1 / (1 + abs(i - j))
      }
    }
  }

  if (type=="triangular") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i,j] <- (abs(j - i)<=spanSamples)*((spanSamples + 1) - abs(j - i))/(spanSamples + 1)
      }
    }
  }

  if (type=="block") {
    for (i in 1:n) {
      for (j in 1:n) {
        W[i,j] <- (abs(j - i)<=spanSamples)*1
      }
    }
  }

  if (plot) {
    require(gplots)
    suppressWarnings(heatmap.2(W,
              trace="none",
              Rowv=FALSE,
              Colv = FALSE,
              col = colorRampPalette(c("black","yellow")),
              main = paste(type,"weight matrix with span", span)))
  }

  return(W)
}

##############################################

#' The DCARS function performs testing for differential correlation across ranked samples.
#'
#' @title DCARS
#' @param dat a genes x samples gene expression rank matrix, should be already converted to ranks with first column lowest survival and last column highest survival
#' @param xname name of row of dat to test together with yname
#' @param yname name of row of dat to test together with xname
#' @param W weight matrix for weighted correlations,
#' @param rangeMin minimum range of weighted correlation vector to include for permutation testing
#' @param wcormin minimum absolute value weighted correlation vector to include for permutation testing
#' @param statmin minimum value DCARS test statistic to include for permutation testing
#' @param weightedConcordanceFunction concordance function to use, defaults to weighted Pearson correlation. User can provide their own function, with arguments fun(x,y,w).
#' @param weightedConcordanceFunctionW either "vector" or "matrix", determining if the function given in weightedConcordanceFunction argument takes in a single vector as input for w, or a matrix.
#' @param plot if TRUE plot observed weighted correlatin vector
#' @param niter number of iterations for permutation testing
#' @param extractTestStatisticOnly if TRUE, extract only the DCARS test statistic without permutation testing
#' @param extractWcorSequence if TRUE, extract only the weighted correlation vector without permutation testing
#' @param extractPermutationTestStatistics if TRUE, extract only the DCARS test statistic without permutation testing
#' @param verbose if TRUE, print updates
#' @param ... additional arguments passing on to weightMatrix()

#' @return either single value (p-value or test statistic), vector (local weighted correlation), or list (combination of above) depending on the input parameters.

#' @examples
#'
#' data(STRING)
#' data(SKCM)
#' SKCM_rank = t(apply(SKCM,1,rank))
#'
#' # highly significantly DCARS gene pair: SKP1 and SKP2
#' # calculates p-value based on permutation
#' DCARS(SKCM_rank,"SKP1","SKP2",plot=TRUE)
#' # extract only the test statistic
#' DCARS(SKCM_rank,"SKP1","SKP2", extractTestStatisticOnly = TRUE)
#'
#' # not significantly DCARS gene pair: EIF3C and EIF5B
#' # calculates p-value based on permutation
#' DCARS(SKCM_rank,"EIF3C","EIF5B",plot=TRUE)
#' # extract only the test statistic
#' DCARS(SKCM_rank,"EIF3C","EIF5B", extractTestStatisticOnly = TRUE)
#'
#' # build weight matrix
#' W = weightMatrix(ncol(SKCM_rank), type = "triangular", span = 0.5, plot = TRUE)
#'
#' # extract DCARS test statistics
#' SKCM_stats = DCARSacrossNetwork(SKCM_rank,edgelist = STRING,
#'                                 W = W, extractTestStatisticOnly = TRUE,
#'                                 verbose = FALSE)
#' sort(SKCM_stats,decreasing=TRUE)[1:10]
#'
#' @export

DCARS = function(dat, xname, yname, W = NULL, rangeMin = 0, wcormin = 0,
                 statmin = 0,
                 weightedConcordanceFunction = weightedPearson,
                 weightedConcordanceFunctionW = "vector",
                 extractTestStatisticOnly = FALSE, extractWcorSequenceOnly = FALSE,
                 plot = FALSE, niter = 100, extractPermutationTestStatistics = FALSE,
                 verbose = FALSE, ...)
{
  weightedcor = weightedConcordanceFunction

  if (any(!c(xname, yname) %in% rownames(dat))) {
    message("xname and/or yname are not found in rownames(dat), returning NA")
    return(NA)
  }
  if (verbose) {
    print(paste(xname, yname))
  }
  if (is.null(W)) {
    message("weight matrix not specified, generating weight matrix now")
    W = weightMatrix(ncol(dat), ...)
  }
  if (!weightedConcordanceFunctionW %in% c("vector","matrix")) {
    stop("Please specify either \"vector\" or \"matrix\" for the weightedConcordanceFunction")
  }
  x = dat[xname, ]
  y = dat[yname, ]

  if (weightedConcordanceFunctionW == "vector") {
    wcor = sapply(1:nrow(W), function(i) weightedcor(x, y, W[i, ]))
  } else {
    wcor = weightedcor(x, y, W)
  }

  if (any(!is.finite(wcor))) {
    i = which(!is.finite(wcor))
    wcor[i] <- mean(c(wcor[c(i - 1, i + 1)]))
  }
  if (plot) {
    plot(wcor, type = "l", ylim = c(-1, 1), lwd = 3, col = "blue")
    abline(h = 0, lty = 2)
  }
  stat = sd(wcor)
  if (extractWcorSequenceOnly & extractTestStatisticOnly) {
    message("both extractWcorSequenceOnly and extractTestStatisticOnly have been specified as TRUE, returning both as a list")
    return(list(wcorSequence = wcor, TestStatistic = stat))
  }
  if (extractWcorSequenceOnly)
    return(wcor)
  if (extractTestStatisticOnly)
    return(stat)
  if (stat < statmin)
    return(NA)
  if (diff(range(wcor)) < rangeMin)
    return(NA)
  if (max(abs(wcor)) < wcormin)
    return(NA)
  if (verbose) {
    message("performing permutation test")
  }
  sds = replicate(niter, {
    o = sample(1:length(x))
    xo = x[o]
    yo = y[o]

    if (weightedConcordanceFunctionW == "vector") {
      wcor = sapply(1:nrow(W), function(i) weightedcor(xo, yo, W[i, ]))
    } else {
      wcor = weightedcor(xo, yo, W)
    }

    if (any(!is.finite(wcor))) {
      i = which(!is.finite(wcor))
      wcor[i] <- mean(c(wcor[c(i - 1, i + 1)]))
    }
    print("done one permutation")
    sd(wcor)
  })
  if (extractPermutationTestStatistics) {
    if (verbose) {
      print("extract permuted test statistics")
    }
    return(list(PermutedTestStatistics = sds))
  }
  return(mean(sds >= stat))
}

##############################################

#' performs DCARS method across all edges listed in the network for the (already ranked) genes x samples matrix dat
#'
#' @title DCARSacrossNetwork
#' @param dat a genes x samples gene expression rank matrix, should be already converted to ranks with first column lowest survival and last column highest survival
#' @param edgelist is a 2 column character matrix with the genes to test for DCARS. if edgelist has more than 2 columns then the first two are taken
#' @param edgeNames is the name to assign to each edge, defaults to rownames(edgelist). if edgelist has no rownames then defaults to two gene names pasted with "_" in alphabetical order
#' @param ... additional parameters passed to DCARS()
#' @return value of this function depends on arguments passed into the DCARS() function (e.g. if extractTestStatisticOnly is set as TRUE)

#' @examples
#' data(STRING)
#' data(SKCM)
#' SKCM_rank = t(apply(SKCM,1,rank))
#'
#' # highly significantly DCARS gene pair: SKP1 and SKP2
#' # calculates p-value based on permutation
#' DCARS(SKCM_rank,"SKP1","SKP2",plot=TRUE)
#' # extract only the test statistic
#' DCARS(SKCM_rank,"SKP1","SKP2", extractTestStatisticOnly = TRUE)
#'
#' # not significantly DCARS gene pair: EIF3C and EIF5B
#' # calculates p-value based on permutation
#' DCARS(SKCM_rank,"EIF3C","EIF5B",plot=TRUE)
#' # extract only the test statistic
#' DCARS(SKCM_rank,"EIF3C","EIF5B", extractTestStatisticOnly = TRUE)
#'
#' # build weight matrix
#' W = weightMatrix(ncol(SKCM_rank), type = "triangular", span = 0.5, plot = TRUE)
#'
#' # extract DCARS test statistics
#' SKCM_stats = DCARSacrossNetwork(SKCM_rank,edgelist = STRING,
#'                                 W = W, extractTestStatisticOnly = TRUE,
#'                                 verbose = FALSE)
#' sort(SKCM_stats,decreasing=TRUE)[1:10]
#' @export

DCARSacrossNetwork = function(dat,edgelist,edgeNames = rownames(edgelist),...) {
  # performs DCARS method across all edges listed in the network
  # for the (already ranked) genes x samples matrix dat
  # value of this function depends on arguments passed into the DCARS() function
  # (e.g. if extractTestStatisticOnly is set as TRUE)
  # dat: is a genes x samples gene expression rank matrix, should be already converted to ranks
  # with first column lowest survival and last column highest survival
  # edgelist: is a 2 column character matrix with the genes to test for DCARS
  # if edgelist has more than 2 columns then the first two are taken
  # edgeNames: is the name to assign to each edge, defaults to rownames(edgelist)
  # if edgelist has no rownames then defaults to two gene names pasted with "_" in alphabetical order

  if (is.null(edgeNames)) {
    if (is.null(rownames(edgelist))) {
      edgeNames = apply(edgelist, 1, function(x) paste0(sort(x[1:2]),
                                                        collapse = "_"))
    }
    else {
      edgeNames = rownames(edgelist)
    }
  }
  DCARSresult = apply(edgelist[, 1:2, drop = FALSE], 1, function(x) DCARS(dat = dat,
                                                                          xname = x[1], yname = x[2], ...))
  if (is.matrix(DCARSresult)) {
    colnames(DCARSresult) <- edgeNames
  }
  else {
    names(DCARSresult) <- edgeNames
  }
  return(DCARSresult)
}

##############################################

#' The pospos function calculates the number of observations that are positive in both x and y. Used for calculating weighted Kendall tau measure of association
#'
#' @title pospos
#' @param x numeric vector of non-negative data x
#' @param y numeric vector of non-negative data y
#' @param offset should 1 be added when pospos is 0 or all
#' @return \code{numeric} of weighted correlations for the sequence of weights given

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' pospos(x,y,offset = TRUE)
#'
#' @export

pospos = function(x,y, offset = TRUE) {
  # offset means either adding or subtracting 1 so that the proportion is not 0 or 1.

  pospos = x > 0 & y > 0
  if (offset) {
    if (sum(pospos) == length(pospos)) pospos[1] <- FALSE
    if (sum(pospos) == 0) pospos[1] <- TRUE
  }
  return(mean(pospos)) # small pospos means many zeros
}

##############################################

#' The logistic function calculates the (inverse) logistic for a given p, where p lies between 0 and 1. Large values indicates many zeros.
#'
#' @title logistic
#' @param x numeric vector of non-negative data x
#' @param y numeric vector of non-negative data y
#' @param offset should 1 be added when pospos is 0 or all
#' @return \code{numeric} of weighted correlations for the sequence of weights given

#' @examples
#'
#' p = seq(0.1,0.9, length.out = 100)
#' l = logistic(p)
#' plot(p,l)
#'
#' @export

logistic = function(p) log((1-p)/(p))


##############################################

#' Performs Fisher's Z Transformation test for differential correlation.
#'
#' @title fisherZtransformTest
#' @param dat dat is a gene expression matrix (rows genes columns samples)
#' @param xname the name of the first gene
#' @param yname name of the second gene
#' @param classLabels a vector of two class labels, if NULL then split into 50/50. odd number of samples will ignore the middle sample
#' @return single value for p-value

#' @examples
#'
#' data(STRING)
#' data(SKCM)
#' fisherZtransformTest(dat = SKCM,xname = "EIF3C",yname = "EIF5B")
#'
#' @export

fisherZtransformTest = function(dat,xname,yname,classLabels=NULL) {
  # Fisher Z-transformation test
  # dat is a gene expression matrix
  # xname is the name of the first gene
  # yname is name of the second gene
  # classLabels is a vector of two class labels, if NULL then split into 50/50
  # odd number of samples will ignore the middle sample

  if (any(!c(xname,yname) %in% rownames(dat))) {
    message("xname and/or yname are not found in rownames(dat), returning NA")
    return(NA)
  }

  classLabelsGiven = !is.null(classLabels)

  if (!classLabelsGiven) {

    classLabels = rep("None",ncol(dat))
    nsamp = floor(ncol(dat)/2)
    classLabels[1:nsamp] <- "Low"
    classLabels[length(classLabels):(length(classLabels)-nsamp+1)] <- "High"
    classes = c("Low","High")

  } else {

    if (length(classLabels)!=ncol(dat)) {
      warning("length of classLabels is not equal to number of columns in dat, returning NA")
      return(NA)
    }

    if (length(unique(classLabels))==1) {
      warning("there is only one class label present, two are needed, returning NA")
      return(NA)
    }

    if (length(unique(classLabels))!=2) {
      warning("there are more than 2 class labels present, taking only the first two of unique class labels")
      classes = unique(classLabels)[1:2]
      print(classes)
    }

    classes = unique(classLabels)

  }

  x = dat[xname,]
  y = dat[yname,]

  x1 = x[classLabels==classes[1]]
  y1 = y[classLabels==classes[1]]
  x2 = x[classLabels==classes[2]]
  y2 = y[classLabels==classes[2]]

  n1 = length(x1)
  n2 = length(x2)

  rc1 = cor(x1,y1)
  rc2 = cor(x2,y2)

  z1 = 0.5*log((1+rc1)/(1-rc1))
  z2 = 0.5*log((1+rc2)/(1-rc2))

  SEdiff = sqrt(1/(n1-3) + 1/(n2-3))

  dz = abs(z1-z2)/SEdiff
  return(2*(1-pnorm(dz)))
}

##############################################

#' Fits a linear model with interaction term and tests for significance of the interaction term.
#'
#' @title LinearModelInteractionTest
#' @param dat dat is a gene expression matrix (rows genes columns samples)
#' @param xname the name of the first gene
#' @param yname name of the second gene
#' @param response a vector of two class labels, if NULL then split into 50/50. odd number of samples will ignore the middle sample
#' @return single value for p-value

#' @examples
#'
#' data(STRING)
#' data(SKCM)
#' LinearModelInteractionTest(dat = SKCM,xname = "EIF3C",yname = "EIF5B")
#'
#' @export

LinearModelInteractionTest = function(dat,xname,yname, response = NULL) {
  # Fitting linear model and testing for interaction effect
  # dat is a gene expression matrix. rows are genes and columns samples
  # xname is the name of the first gene
  # yname is name of the second gene
  # response = continuous response, e.g. survival information (assuming all are uncensored)
  # returns the p-value from this test

  if (any(!c(xname,yname) %in% rownames(dat))) {
    message("xname and/or yname are not found in rownames(dat), returning NA")
    return(NA)
  }

  x = dat[xname,]
  y = dat[yname,]

  if (is.null(response)) {
    response = 1:length(x)
  }

  fit = lm(response ~ x*y) # this tests the interaction effect between x and y
  return(summary(fit)$coef[4,4])

}

##############################################

#' Calculates weighted Pearson correlation between x and y.
#' This function has been superceded by weightedPearson() and
#' similar matrix function weightedPearson_matrix()
#'
#' @title weightedcor
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @param w weight vector
#' @return weighted correlation value between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' w = runif(100)
#' weightedcor(x,y,w)
#'
#' @export


weightedcor = function(x,y,w) {

  # calculates weighted Pearson correlation between x and y
  # x and y are data vectors
  # w is a weight vector

  nw = sum(w)
  wssx = nw*sum(w*(x^2)) - sum(w*x)^2
  wssy = nw*sum(w*(y^2)) - sum(w*y)^2
  wssxy = nw*sum(w*x*y) - sum(w*x)*sum(w*y)
  wcor = wssxy/sqrt(wssx*wssy)
  return(wcor)
}

##############################################

#' the weightedPearson function
#'
#' @title weightedPearson
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @param w weight vector, values should be between 0 and 1
#' @return \code{numeric} weighted correlation value between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' w = runif(100)
#' weightedPearson(x,y,w)
#'
#' @export

weightedPearson = function(x, y, w) {
  nw = sum(w)
  wssx = nw * sum(w * (x^2)) - sum(w * x)^2
  wssy = nw * sum(w * (y^2)) - sum(w * y)^2
  wssxy = nw * sum(w * x * y) - sum(w * x) * sum(w * y)
  wcor = wssxy/sqrt(wssx * wssy)
  return(wcor)
}

##############################################

#' The weightedPearson_matrix function calculates a vector of weighted correlations for two given data vectors, for a matrix of given weights.
#'
#' @title weightedPearson_matrix
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @param W weight matrix, values should be between 0 and 1, number of columns should be the same as length(x) and length(y)
#' @return \code{vector} weighted correlation values between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' W = weightMatrix(100)
#' weightedPearson_matrix(x,y,w)
#'
#' @export

weightedPearson_matrix = function(x, y, W) {
  x2 = x^2
  y2 = y^2
  xy = x*y

  wcorVec = rep(0,nrow(W))

  for (i in 1:nrow(W)) {

    w = W[i,]

    wx = w*x
    wy = w*y

    nw = sum(w)
    wssx = nw * sum(w * x2) - sum(wx)^2
    wssy = nw * sum(w * y2) - sum(wy)^2
    wssxy = nw * sum(w * xy) - sum(wx) * sum(wy)
    wcor = wssxy/sqrt(wssx * wssy)

    wcorVec[i] <- wcor
  }
  return(wcorVec)
}

##############################################

#' the weightedSpearman function
#'
#' @title weightedSpearman
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @param w weight vector, values should be between 0 and 1
#' @return \code{numeric} weighted correlation value between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' w = runif(100)
#' weightedSpearman(x,y,w)
#'
#' @export

weightedSpearman = function(x,y,w) {
  xr = rank(x)
  yr = rank(y)
  return(weightedPearson(xr,yr,w))
}

##############################################

#' The weightedSpearman_matrix function calculates a vector of weighted correlations for two given data vectors, for a matrix of given weights.
#'
#' @title weightedSpearman_matrix
#' @param x x and y are data vectors
#' @param y x and y are data vectors
#' @param W weight matrix, values should be between 0 and 1, number of columns should be the same as length(x) and length(y)
#' @return \code{vector} weighted correlation values between x and y

#' @examples
#'
#' x = rnorm(100)
#' y = rnorm(100)
#' W = weightMatrix(100)
#' weightedSpearman_matrix(x,y,w)
#'
#' @export

weightedSpearman_matrix = function(x,y,W) {
  xr = rank(x)
  yr = rank(y)
  return(weightedPearson_matrix(xr,yr,W))
}

##############################################

#' the weightedKendallStar function calculates weighted Tau*, where Tau* is described in Pimentel et al (2015) doi:10.1016/j.spl.2014.09.002. This association measure is defined for zero-inflated, non-negative random variables.
#'
#' @title weightedKendallStar
#' @param x x and y are non-negative data vectors
#' @param y x and y are non-negative data vectors
#' @param w weight vector, values should be between 0 and 1
#' @return \code{numeric} weighted Tau* association value between x and y

#' @examples
#'
#' x = pmin(0,rnorm(100))
#' y = pmin(0,rnorm(100))
#' w = runif(100)
#' weightedKendallStar(x,y,w)
#'
#' @export

weightedKendallStar = function(x, y, w) {

  if (any(x < 0 | y < 0)) stop("x and/or y values have negative values")

  posx = x > 0
  posy = y > 0

  pospos = posx & posy

  x_diff = outer(x[pospos],x[pospos], FUN = "-")
  y_diff = outer(y[pospos],y[pospos], FUN = "-")
  w_diff = outer(w[pospos],w[pospos])
  w_diff[lower.tri(w_diff, diag = TRUE)] <- 0

  p_conc = sum((x_diff*y_diff*upper.tri(x_diff) > 0)*w_diff)/sum(w_diff)
  p_disc = sum((x_diff*y_diff*upper.tri(x_diff) < 0)*w_diff)/sum(w_diff)

  tau_11 = p_conc - p_disc

  p_11 = sum(w*pospos)/sum(w)

  p_00 = sum(w*(!posx & !posy))/sum(w)

  p_01 = sum(w*(!posx & posy))/sum(w)

  p_10 = sum(w*(posx & !posy))/sum(w)

  # define p_1 and p_2
  # p_1 (weighted) proportion of times the x-values with y = 0 is greater than
  # the x-values with y > 0
  x_10 = x[!posy]
  x_11 = x[posy]

  w_x_10 = w[!posy]
  w_x_11 = w[posy]

  w_x_outer = outer(w_x_10, w_x_11)

  if (sum(w_x_outer) == 0) {
    p_1 <- 0
  } else {
    p_1 = sum(outer(x_10, x_11, FUN = ">")*w_x_outer) / sum(w_x_outer)
  }

  # analogous for y
  y_10 = y[!posx]
  y_11 = y[posx]

  w_y_10 = w[!posx]
  w_y_11 = w[posx]

  w_y_outer = outer(w_y_10, w_y_11)

  if (sum(w_y_outer) == 0) {
    p_2 <- 0
  } else {
    p_2 = sum(outer(y_10, y_11, FUN = ">")*w_y_outer) / sum(w_y_outer)
  }

  tauStar = (p_11^2)*tau_11 +
    2*(p_00*p_11 - p_01*p_10) +
    2*p_11*(p_10*(1-2*p_1) + p_01*(1-2*p_2))

  return(tauStar)
}

##############################################

#' the weightedVariance function
#'
#' @title weightedVariance
#' @param x x is a data vector
#' @param y default to NULL, if given it is ignored
#' @param w weight vector, values should be between 0 and 1
#' @return \code{numeric} weighted variance value for x

#' @examples
#'
#' x = rnorm(100)
#' w = runif(100)
#' weightedVariance(x,w)
#'
#' @export

weightedVariance = function(x, y = NULL, w) {
  # args x,w (if y given, it is ignored)

  nw = sum(w)
  wssx = nw * sum(w * (x^2)) - sum(w * x)^2
  return(wssx)

}

##############################################

#' The weightedVariance_matrix function calculates a vector of weighted variances for a given data vector, for a matrix of given weights.
#'
#' @title weightedVariance_matrix
#' @param x x is a data vector
#' @param y default to NULL, if given, it is ignored
#' @param W weight matrix, values should be between 0 and 1, number of columns should be the same as length(x) and length(y)
#' @return \code{vector} weighted variances of x for the weights

#' @examples
#'
#' x = rnorm(100)
#' W = weightMatrix(100)
#' weightedVariance_matrix(x,y,w)
#'
#' @export

weightedVariance_matrix = function(x, y = NULL, W) {
  # args x,w (if y given, it is ignored)

  x2 = x^2

  wssxVec = rep(0,nrow(W))

  for (i in 1:nrow(W)) {

    w = W[i,]

    nw = sum(w)
    wssx = nw * sum(w * (x2)) - sum(w * x)^2
    wssxVec[i] <- wssx

  }

  return(wssxVec)

}

##############################################

#' The stretch function calculates a 'stretched'/standardised vector of weighted correlations, using previously calculated upper and lower bounds of the statistic (same length vectors). Values are stretched higher or lower symmetrically around 0.
#'
#' @title stretch
#' @param wcor weighted correlation vector
#' @param upper either single numeric or vector of same length as wcor for upper bound of correlation. Default 1
#' @param lower either single numeric or vector of same length as wcor for lower bound of correlation. Default -1
#' @return \code{vector} standardised weighted correlation vector

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' w = runif(100)
#' wcor = weightedKendallStar(x,y,w)
#' stretch(wcor, upper = 0.8, lower = -0.3)
#'
#' @export

stretch = function(wcor,upper = 1,lower = -1) {
  # given a vector of observed weighted correlations,
  # and previously calculated upper and lower bounds of the statistic (same length vectors)
  # calculate a 'stretched'/standardised vector of weighted correlations

  if (length(upper) == 1 & length(lower) == 1 & upper[1] == 1 & lower[1] == -1) return(wcor)

  if (length(wcor) != length(upper) | length(wcor) != length(lower)) stop("Need upper and lower bounds to have same length as weighted correlation")

  swcor = rep(0,length(wcor))

  pos = wcor > 0 # logical vector of positive wcors
  neg = wcor < 0 # logical vector of negative wcors

  swcor[pos] <- wcor[pos]/upper[pos] # stretch upwards

  swcor[neg] <- wcor[neg]/(-lower[neg]) # stretch downwards, retaining sign

  return(swcor)
}

##############################################

#' The boundsKendallStar function calculates the upper and lower bounds the weighted zero-inflated Kendall's tau star association measure, given a matrix containing a number of weights.
#'
#' @title boundsKendallStar
#' @param x data vector
#' @param y data vector
#' @param W weight matrix, number of columns correspond to length of x and y.
#' @return \code{list} list with two objects named "upper" and "lower", vectors of upper and lower bounds on the association measure

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' w = runif(100)
#' wcor = weightedKendallStar(x,y,w)
#' boundsKendallStar(x,y,w)
#'
#' @export

boundsKendallStar = function(x,y,W) {
  # for data vectors x and y, calculate the upper and lower bounds
  # on the weighted zero-inflated Kendall's tau star association measure
  # given a matrix containing a number of weights

  # output is a list containing vectors of upper and lower bounds
  # on the association measure

  # note this only depends on the weighted proportion of
  # samples with zero values

  require(Matrix)

  if (!is.matrix(W)) {
    W <- t(as.matrix(W))
  }

  x_zero = 1*(x == 0)
  y_zero = 1*(y == 0)
  rowSums_W = Matrix::rowSums(W)

  p_x_weighted = ( W %*% x_zero ) / rowSums_W
  p_y_weighted = ( W %*% y_zero ) / rowSums_W

  bound_upper = 1 - pmax(p_x_weighted^2, p_y_weighted^2)
  bound_lower = pmax(0, (1 - p_x_weighted - p_y_weighted))^2 - 2*(1-p_x_weighted)*(1-p_y_weighted)

  return(list(lower = c(bound_lower),
              upper = c(bound_upper)
  ))

}

##############################################

#' The sweightedKendallStar function is a convenience function to get swcor for weighted zero-inflated kendall's tau
#'
#' @title sweightedKendallStar
#' @param x data vector
#' @param y data vector
#' @param w weights vector
#' @return \code{vector} of stretched weighted correlations using zero-inflated kendall's tau association measure

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' w = runif(100)
#' wcor = weightedKendallStar(x,y,w)
#' bounds = boundsKendallStar(x,y,w)
#' stretch(wcor, upper = bounds[["upper"]], lower = bounds[["lower"]]
#'
#' sweightedKendallStar(x,y,w)
#'
#' @export

sweightedKendallStar = function(x,y,w) {
  # convenience function to get swcor for weighted zero-inflated kendall's tau

  wcor = weightedKendallStar(x,y,w)

  bounds = boundsKendallStar(x,y,w)

  swcor = stretch(wcor,upper = bounds[["upper"]], lower = bounds[["lower"]])

  return(swcor)
}

##############################################

#' The thin function extracts the rows of a matrix evenly so that roughly n number of rows remain. Used for thinning down the weight matrix to speed up overall computation.
#'
#' @title thin
#' @param W matrix
#' @param n rough number of rows to keep
#' @return \code{matrix} of thinned matrix keeping only roughly n rows.

#' @examples
#'
#' W = weightMatrix(500)
#' W_small = thin(W, n = 100)
#'
#' @export

thin = function(W, n = 100) {
  # given a matrix W, extract the rows so that
  # the total number of rows is roughly n

  N = nrow(W)

  if (n == N) return(W)

  if (n > N) stop("need to set n to less than nrow(W)")

  new_n = floor(N/n)

  index = seq(from = 1, to = N, by = new_n)

  return(W[index,])
}

##############################################

#' The weightedStretchedKendallStar_matrix function calculates the stretched weighted Tau*, where where Tau* is described in Pimentel et al (2015) doi: 10.1016/j.spl.2014.09.002, and upper and lower bounds given in Denuit et al (2017) doi: 10.1016/j.spl.2017.03.005
#'
#' @title weightedStretchedKendallStar_matrix
#' @param x numeric vector of non-negative data x
#' @param y numeric vector of non-negative data y
#' @param W weight matrix, same columns as length x and y
#' @param stretch logical default TRUE, whether to stretch based on calculated upper and lower bounds
#' @return \code{vector} of weighted correlations for the sequence of weights given

#' @examples
#'
#' x = pmax(0,rnorm(100))
#' y = pmax(0,rnorm(100))
#' W = weightMatrix(100)
#' wcor = weightedStretchedKendallStar_matrix(x, y, W)
#'
#' @export

weightedStretchedKendallStar_matrix = function(x, y, W, stretch = TRUE) {

  # function to calculate (stretched) weighted Tau*, where Tau* is described in
  # Pimentel et al (2015)
  # this association measure is defined on zero-inflated, non-negative random variables.
  # W is a weight matrix with columns equal to the number of samples

  if (any(x < 0 | y < 0)) stop("x and/or y values have negative values")

  posx = x > 0
  posy = y > 0

  pospos = posx & posy

  # look into calculating weighted tau from the positive values
  # tau_11 = cor(x[pospos], y[pospos], method = "kendall")

  x_diff = outer(x[pospos],x[pospos], FUN = "-")
  y_diff = outer(y[pospos],y[pospos], FUN = "-")
  xy_diff_upper_tri_x_diff = x_diff*y_diff*upper.tri(x_diff)

  x_10 = x[!posy]
  x_11 = x[posy]

  y_10 = y[!posx]
  y_11 = y[posx]

  # weight specific quantities
  tauStarVec = rep(0,nrow(W))
  for (i in 1:nrow(W)) {

    w = W[i,]

    w_diff = outer(w[pospos],w[pospos])
    w_diff[lower.tri(w_diff, diag = TRUE)] <- 0

    sum_w = sum(w)

    sum_w_diff = sum(w_diff)
    if (sum_w_diff != 0) {
      p_conc = sum((xy_diff_upper_tri_x_diff > 0)*w_diff)/sum_w_diff
      p_disc = sum((xy_diff_upper_tri_x_diff < 0)*w_diff)/sum_w_diff

      tau_11 = p_conc - p_disc
    } else {
      tau_11 = 0
    }

    p_11 = sum(w*pospos)/sum_w

    p_00 = sum(w*(!posx & !posy))/sum_w

    p_01 = sum(w*(!posx & posy))/sum_w

    p_10 = sum(w*(posx & !posy))/sum_w

    # define p_1 and p_2
    # p_1 (weighted) proportion of times the x-values with y = 0 is greater than
    # the x-values with y > 0
    w_x_10 = w[!posy]
    w_x_11 = w[posy]

    w_x_outer = outer(w_x_10, w_x_11)

    if (sum(w_x_outer) == 0) {
      p_1 <- 0
    } else {
      p_1 = sum(outer(x_10, x_11, FUN = ">")*w_x_outer) / sum(w_x_outer)
    }

    # analogous for y
    w_y_10 = w[!posx]
    w_y_11 = w[posx]

    w_y_outer = outer(w_y_10, w_y_11)

    if (sum(w_y_outer) == 0) {
      p_2 <- 0
    } else {
      p_2 = sum(outer(y_10, y_11, FUN = ">")*w_y_outer) / sum(w_y_outer)
    }

    tauStar = (p_11^2)*tau_11 +
      2*(p_00*p_11 - p_01*p_10) +
      2*p_11*(p_10*(1-2*p_1) + p_01*(1-2*p_2))

    tauStarVec[i] <- tauStar
  }

  if (stretch) {
    bounds = boundsKendallStar(x,y,W)
    tauStarVec = stretch(tauStarVec,upper = bounds[["upper"]], lower = bounds[["lower"]])
  }

  return(tauStarVec)
}

##############################################

#' the stratifiedSample function takes in a vector of statistics, and returns a vector of indices for the statistics that should be used for permutation testing
#'
#' @title stratifiedSample
#' @param stats a vector of statistics
#' @param length rough number of indices to extract
#' @return \code{vector} a vector of indices for the statistics that should be used for permutation testing

#' @examples
#'
#'stats = rnorm(5000)
#'sampleindices = stratifiedSample(stats)
#'length(sampleindices)
#'
#' @export
#'
stratifiedSample = function(stats, length = 100) {
  nsamps = 3
  ranges = range(stats[is.finite(stats)])
  fac = cut(stats, seq(from = ranges[1] - 1e-4, to = ranges[2], length.out = ceiling(length/nsamps)))
  sampleindices = unlist(tapply(1:length(stats), fac, function(x) {
    if (length(x) < nsamps)
      return(x)
    sample(x, nsamps, replace = FALSE)
  }))
  sampleindices = unique(sampleindices[!is.na(sampleindices)])
  return(sampleindices)
}

##############################################

#' the getLoessCriticalValue function
#'
#' @title getLoessCriticalValue
#' @param stats a vector of DCARS test statistics for which permutation testing has been done
#' @param pvals a vector of pvals for DCARS test
#' @param signifValue significance value (default 0.05)
#' @param plot logical (default FALSE)
#' @return \code{numeric} a critical value for the unadjusted significance value

#' @examples
#'
#'stats = rchisq(100,1)
#'pvals = 1 - pchisq(stats,1)
#'criticalValue = getLoessCriticalValue(stats, pvals, plot = TRUE)
#'criticalValue
#'
#' @export
#'
getLoessCriticalValue = function(stats, pvals, signifValue = 0.05, plot = FALSE) {
  if (length(stats) != length(pvals)) stop("stats and pvals vectors should be equal in length!")

  lw = loess(pvals~stats)
  # points(x[1:ceiling(0.9*nrow(x)),1],lw$fitted[1:ceiling(0.9*nrow(x))], type = "l", col = "red")
  lwval = which(sapply(1:(length(stats)-1),function(i)lw$fitted[i]>signifValue & lw$fitted[i+1]<signifValue))[1]
  criticalValue = as.numeric(stats[lwval])
  if (plot) {
    plot(stats, pvals)
    points(stats[1:ceiling(0.9*length(stats))],lw$fitted[1:ceiling(0.9*length(stats))], type = "l", col = "red")
    abline(v = criticalValue, h = signifValue, col = "blue", lwd = 2, lty = 2)
    legend("topright", bty = "n", legend = paste0("Critical value for P = ",signifValue, ": ", signif(criticalValue, 3)))
  }
  return(criticalValue)
}

##############################################

#' the plotColouredExpression function plots a 3 panel scatterplot of the gene pairs split by early, mid, and late in the sample ordering.
#'
#' @title plotColouredExpression
#' @param branchData is a list containing matrices of the cell expression per branch, assumed that the columns of each matrix in branchData is ordered by pseudotime. If branchData is not given as a list, it will be converted into a list containing branchData.
#' @param genepair is either a single character string with an underscore, or a two length character vector
#' @param subsetBranch subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
#' @return \code{ggplot} a ggplot object of scatterplots of expression split by sample ordering

#' @examples
#'
#'
#'
#'
#'
#'
#' @export
#'
plotColouredExpression = function(branchData, genepair, subsetBranch = NULL, n = 3, fittedline = TRUE) {
  require(ggplot2)
  if (length(genepair) == 1) {
    genepair = unlist(strsplit(genepair, "_"))
  }
  else {
    genepair = genepair[1:2]
  }
  if (!is.list(branchData)) {
    branchData = list(Branch = branchData)
  }
  if (is.null(names(branchData))) {
    names(branchData) <- paste0("Branch_", 1:length(branchData))
  }
  gdf = do.call(rbind, lapply(branchData, function(branch) {
    gdf_list_1 = data.frame(Sample = colnames(branch), order = 1:ncol(branch),
                            ExpressionGene1 = branch[genepair[1], ], ExpressionGene2 = branch[genepair[2],
                                                                                              ])
    if (n > 1) {
      gdf_list_1$ordercut = cut(gdf_list_1$order, n, labels = unlist(ifelse(n == 3, list(c("Early",
                                                                                           "Middle", "Late")), list(paste0("Group ", 1:n)))))
    } else {
      gdf_list_1$ordercut = rep("All cells", times = nrow(gdf_list_1))
    }
    return(gdf_list_1)
  }))
  gdf$branch = rep(names(branchData), times = unlist(lapply(branchData,
                                                            ncol)))
  if (!is.null(subsetBranch)) {
    gdf_sub = subset(gdf, branch %in% subsetBranch)
    if (nrow(gdf_sub) == 0)
      stop("no branches with names in subsetBranch, please re-run with correct names (should match names of branchData)")
  }
  else {
    gdf_sub = gdf
  }
  g = ggplot(gdf_sub, aes(x = ExpressionGene1, y = ExpressionGene2,
                          color = order)) +
    geom_point(show.legend = FALSE, alpha = 0.7) +
    facet_grid(branch ~ ordercut, scales = "free_y") +
    scale_color_gradientn(colours = c("orange", "blue")) +
    xlab(genepair[1]) + ylab(genepair[2]) + theme_minimal() +
    NULL
  if (fittedline) {
    g = g +
      geom_smooth(colour = "black", fill = NA, linetype = "dashed",
                  method = "lm") +
      NULL
  }
  return(g)
}

##############################################

#' the plotEgoNetwork function
#'
#' @title plotEgoNetwork plots network graphs with edges coloured by weights in the network
#' @param hubnode is a character vector of node(s) to include as hub nodes
#' @param g is an igraph network, with E(g)[[weight]] given as DCARS test statistics
#' @param subset is a logical asking if you should subset based on the weight (default FALSE)
#' @param thresh is the subset weight threshold
#' @return \code{igraph} object containing the network graphed. Produces an igraph plot

#' @examples
#'
#'
#'
#'
#'
#'
#' @export
#'
plotEgoNetwork = function(hubnode, network, weight = "weight", subset = FALSE, thresh = NULL) {

  require(igraph)

  # hubnode is a character vector of node(s) to include as hub nodes
  # g is an igraph network, with E(g)[[weight]] given as DCARS test statistics
  # weight is a character vector containing edge weights, associated with different branches
  # subset is a logical asking if you should subset based on the weight
  # thresh is the subset weight threshold

  nodes = unique(names(unlist(neighborhood(network, nodes = hubnode))))
  subego = induced.subgraph(network, vids = nodes)
  subego = simplify(subego, edge.attr.comb="mean")

  if (subset) {
    if (!all(weight %in% names(edge_attr(subego)))) {
      stop("at least one weight missing from edge attributes, either re-specify weights or rerun with subset = FALSE")
    }
    if (is.null(thresh)) {
      message("no threshold given, using 0 as default")
      thresh = 0
    }
    keepedges = apply(sapply(weight, function(w) edge_attr(subego)[[w]] > thresh,
                             simplify = TRUE),1,any)
    subego = subgraph.edges(subego, which(keepedges))
  }

  V(subego)$color = "beige"
  V(subego)$label.color = "black"
  V(subego)$label.family = "sans"
  V(subego)$label.cex = 0.7
  V(subego)$size = 20
  V(subego)$frame.color = "black"
  V(subego)$frame.size = 5
  lyout = layout.davidson.harel(subego)
  width = 5

  maxval = ceiling(max(50*unlist(edge_attr(subego)[weight])))
  colvals = colorRampPalette(c("grey","red"))(maxval)

  par(mfrow=c(1,length(weight)))
  for (i in weight) {
    plot(subego, layout = lyout,
         edge.color = colvals[ceiling(50*edge_attr(subego)[[i]])],
         edge.width = width,
         xlab = i,
         main = hubnode)
  }
  return(subego)
}

##############################################

#' the plotWCorLine function
#'
#' @title plotWCorLine plots weighted correlation vectors as line plots
#' @param wcorsList is a list of matrices, with each matrix gene pair x samples weighted correlation vectors, assumed that they have same number of rows
#' @param gene is either a logical vector matching rows of entries in wcorsList, or a character of a gene
#' @return \code{ggplot} object with line plots

#' @examples
#'
#'
#'
#'
#'
#'
#' @export
#'
plotWCorLine = function(wcorsList, gene) {

  # wcorsList is a list of matrices, with each matrix gene pair x samples weighted correlation vectors, assumed that they have same number of rows
  # gene is either a logical vector matching rows of entries in wcorsList, or a character of a gene
  # matchExact matches gene names by splitting instead of using grep, but is slower

  require(reshape)

  if (class(wcorsList) != "list") {
    wcorsList = list(Branch = wcorsList)
  }

  if (is.null(names(wcorsList))) {
    names(wcorsList) <- paste0("Branch_",1:length(wcorsList))
  }

  if (is.logical(gene[1])) {

    if (length(unique(unlist(lapply(wcorsList,nrow)))) > 1) {
      stop("cannot use logical subset when weighted correlation matrices have differing rows")
    }

    if (length(gene) != nrow(wcorsList[[1]])) {
      stop("cannot use logical subset when length of gene doesn't match nrow of wcorsList matrices")
    }

    wcors_longList = lapply(wcorsList,function(branch){
      melt(t(branch[gene,]))
    })

    gene = ""
  } else {

    gene = paste0(sort(gene), collapse = "|")

    wcors_longList = lapply(wcorsList,function(branch){
      melt(t(branch[grepl(gene,rownames(branch)),]))
    })

  }

  branch_long = do.call(rbind,wcors_longList)
  branch_long = cbind(
    rep(names(wcors_longList), unlist(lapply(wcors_longList,nrow))),
    branch_long)
  colnames(branch_long) = c("branch","SampleOrder", "GenePair","WeightedCorrelation")


  g = ggplot(branch_long,
             aes(x = SampleOrder, y = WeightedCorrelation, group = GenePair, col = GenePair)) +
    geom_line(size = 2, alpha = 0.6) +
    facet_grid(~branch, scales = "free_x") +
    theme_minimal() +
    ylim(c(-1,1)) +
    geom_hline(yintercept = 0, size = 1, colour = "grey") +
    ggtitle(gene) +
    NULL

  return(g)
}

##############################################

#' the plotOrderedExpression function
#'
#' @title plotOrderedExpression plots expression vectors along branches and genes as ribbon plots
#' @param branchData is a list containing matrices of the cell expression per branch, assumed that the columns of each matrix in branchData is ordered by pseudotime. If branchData is not given as a list, it will be converted into a list containing branchData.
#' @param gene is either a single character string with an underscore, or a two length character vector
#' @param xvals is a list containing the x-values associated with the samples in branchData (if NULL, samples will just be plotted against their rank)
#' @param subsetBranch subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
#' @param facet can either be FALSE, "branch", "gene", or "both"
#' @return \code{ggplot} a ggplot object for ribbon plot with points
#'
#' @examples
#'
#'
#'
#'
#' @export
#'
plotOrderedExpression = function(branchData, gene, xvals = NULL, subsetBranch = NULL, facet = FALSE) {

  # branchData is a list containing matrices of the cell expression per branch
  # assumed that the columns of each matrix in branchData is ordered by pseudotime
  # if branchData is not given as a list, it will be converted into a list containing branchData
  # gene is either a single character string with an underscore, or a two length character vector
  # xvals is a list containing the x-values associated with the samples in branchData (if NULL, samples will just be plotted against their rank)
  # subsetBranch is a character vector containing the names of the branches to be plotted. If NULL it will plot all branches
  # facet can either be FALSE, "branch", "gene", or "both"

  require(ggplot2)

  if (!is.list(branchData)) {
    branchData = list(Branch = branchData)
  }

  if (is.null(names(branchData))) {
    names(branchData) <- paste0("Branch_",1:length(branchData))
  }

  gdf_list = sapply(gene, function(g) {

    gdf = do.call(rbind,lapply(branchData, function(branch){
      gdf_list_1 = data.frame(
        Sample = colnames(branch),
        order = 1:ncol(branch),
        ExpressionGene = branch[g,],
        gene = g
      )
      return(gdf_list_1)
    }))
    gdf$branch = rep(names(branchData), times = unlist(lapply(branchData, ncol)))
    return(gdf)
  }, simplify = FALSE)

  gdf = do.call(rbind,gdf_list)

  if (!is.null(subsetBranch)) {
    gdf_sub = subset(gdf, branch %in% subsetBranch)
    if (nrow(gdf_sub) == 0) stop("no branches with names in subsetBranch, please re-run with correct names (should match names of branchData)")
  } else {
    gdf_sub = gdf
  }

  if (!is.null(xvals)) {
    xval = apply(as.matrix(gdf_sub), 1, function(x)xvals[[x["branch"]]][as.numeric(x["order"])])
    gdf_sub$order <- xval
  }

  g = ggplot(gdf_sub, aes(x = order, y = ExpressionGene, colour = gene, fill = gene, linetype = branch, shape = branch)) +
    geom_point() +
    labs(fill = "Gene", col = "Branch") +
    theme_minimal() + geom_smooth() +
    ggtitle(paste0(gene, collapse = ", ")) +
    NULL

  if (facet == "branch") {
    g = g + facet_grid(~branch)
  }

  if (facet == "gene") {
    g = g + facet_grid(gene~.)
  }

  if (facet == "both") {
    g = g + facet_grid(gene~branch)
  }

  return(g)
}


##############################################

#' the plotNetworkPathway function
#'
#' @title plotNetworkPathway plots network graphs with most represented pathways labelled
#' @param sigPairsList list of significant pairs for branches, in igraph network form
#' @param minCommunity minimum number of genes in community to annotate most represented
#' @param pathways is a named list containing gene sets for querying, if left NULL defaults to REACTOME c2.all.v6.1.symbols.gmt downloaded from MSigDB
#' @param geneUniverse is the universe of genes measured, if null its taken as all the sigpairs
#' @return \code{plot} plots of igraph objects
#'
#' @examples
#'
#'
#'
#'
#' @export
#'
plotNetworkPathway = function(sigPairsList, minCommunity = 10, pathways = NULL, geneUniverse = NULL) {

  require(igraph)

  # sigPairsList list of significant pairs for branches, in igraph network form
  # minCommunity minimum number of genes in community to annotate most represented
  # pathways is a named list containing gene sets for querying, if left NULL defaults to REACTOME c2.all.v6.1.symbols.gmt downloaded from MSigDB
  # geneUniverse is the universe of genes measured, if null its taken as all the sigpairs

  if (is.null(pathways)) {
    data(hallmarks)
  } else {
    hallmarks = pathways
  }

  allhallmarks = unique(unlist(hallmarks))

  if (class(sigPairsList) != "list") {
    sigPairsList = list(Branch = sigPairsList)
  }

  if (is.null(geneUniverse)) {
    geneUniverse = unique(unlist(lapply(sigPairsList, function(net) V(net)$name)))
  }

  lyList = lapply(sigPairsList, function(net) {
    ly = layout.auto(net)
    ly = apply(ly,2,function(x) {
      2*(-0.5 + (x - min(x))/(max(x) - min(x)))
    })
    return(ly)
  })

  cmList = lapply(sigPairsList, walktrap.community)
  cmmembershipList = lapply(cmList,membership)
  cm_coordsList = sapply(names(cmList), function(branch) {
    cbind(tapply(lyList[[branch]][,1],cmmembershipList[[branch]],mean),
          tapply(lyList[[branch]][,2],cmmembershipList[[branch]],mean))
  }, simplify = FALSE)
  cm_nList = sapply(names(cmList), function(branch) {
    tapply(lyList[[branch]][,1],cmmembershipList[[branch]],length)
  }, simplify = FALSE)

  message("Calculating most highly represented pathways")

  cm_mostrepresentedList = sapply(names(cmList), function(branch) {

    cm_mostrepresented = sapply(unique(cmmembershipList[[branch]]), function(i){
      # print(i)
      genes = intersect(allhallmarks,
                        toupper(names(cmmembershipList[[branch]][cmmembershipList[[branch]] == i])))
      if (length(genes) == 0) return("")
      allgenes = intersect(allhallmarks, toupper(geneUniverse))

      # calculate the smallest set with the highest membership proportion
      hallmarksn = unlist(lapply(hallmarks, function(x)length(intersect(x,allgenes))))
      hallmarksmembership = unlist(lapply(hallmarks,function(hallmark){
        mean(genes %in% intersect(hallmark,allgenes))
      }))
      if (sum(hallmarksmembership == max(hallmarksmembership)) == 1) {
        return(names(sort(hallmarksmembership, decreasing = TRUE)[1]))
      } else {
        # hallmarksmembershipmax = hallmarksmembership[hallmarksmembership == max(hallmarksmembership)]
        hallmarksnmax = hallmarksn[hallmarksmembership == max(hallmarksmembership)]
        return(names(sort(hallmarksnmax)[1]))
      }
    })
    names(cm_mostrepresented) <- unique(cmmembershipList[[branch]])
    return(cm_mostrepresented)
  }, simplify = FALSE)

  message("Plotting network graphs")

  par(mfrow=c(length(cmList),2))

  for (i in names(cmList)) {

    plot(sigPairsList[[i]],vertex.size = 0, vertex.label.cex = 0.8, vertex.color = "grey",
         edge.width = 2, layout = lyList[[i]], ylab = i,
         vertex.label.color = "black")

    plot(sigPairsList[[i]],vertex.size = 2, vertex.label.cex = 0.01, vertex.color = "grey",
         edge.width = 2, layout = lyList[[i]],
         vertex.label.color = "black")

    text(cm_coordsList[[i]][cm_nList[[i]] >= minCommunity,1],
         cm_coordsList[[i]][cm_nList[[i]] >= minCommunity,2],
         cm_mostrepresentedList[[i]][cm_nList[[i]] >= minCommunity],
         cex = 0.7)
  }

}

##############################################

#' the plotWcorsClusterPathway function
#'
#' @title plotWcorsClusterPathway plots network graphs with most represented pathways labelled
#' @param sigWcorsList named list of weighted correlation matrices per branch, rows are gene pairs and columns are ordered samples, if it is not a list it will be converted into a list
#' @param pathways is a named list containing gene sets for querying, if left NULL defaults to REACTOME c2.all.v6.1.symbols.gmt downloaded from MSigDB
#' @param geneUniverse is the universe of genes measured, if null its taken as all the genes in the sigWcorsList
#' @param cutk number of weighted correlation clusters to cut hclust
#' @param topPathways is the top number of pathways to plot
#' @param cluster is a logical if pathway and weighted correlation clustering should be performed
#' @param label_wrap_cut number of letters to wrap for labelling heatmap
#' @return \code{plot} plots of igraph objects
#'
#' @examples
#'
#'
#'
#'
#' @export
#'
plotWcorsClusterPathway = function(sigWcorsList,pathways = hallmarks,geneUniverse = NULL,cutk = 15,topPathways = 2,cluster = FALSE,label_wrap_cut = 20) {

  # sigWcorsList named list of weighted correlation matrices per branch, rows are gene pairs and columns are ordered samples
  # pathways named list of pathways containing gene names
  # geneUniverse genes in the universe for pathway analysis
  # cutk number of weighted correlation clusters to cut hclust
  # topPathways is the top number of pathways to plot
  # cluster = FALSE
  # label_wrap_cut number of letters to wrap

  require(patchwork)
  require(ggplot2)
  require(reshape)

  if (is.null(pathways)) {
    data(hallmarks)
  } else {
    hallmarks = pathways
  }

  allhallmarks = unique(unlist(hallmarks))

  if (class(sigWcorsList) != "list") {
    sigWcorsList = list(Branch = sigWcorsList)
  }

  if (is.null(geneUniverse)) {
    geneUniverse = toupper(unique(unlist(lapply(sigWcorsList, function(wcors) strsplit(rownames(wcors),"_")))))
  }

  # make list of significant wcors matrix

  sig_wcors_hc = lapply(sigWcorsList,function(wcors) {
    cutree(hclust(dist(wcors, method = "euclidean")),cutk)
  })

  message("Extracted weighted correlation clusters")

  # extract genes from the cutk clusters
  clusterGenes = lapply(sig_wcors_hc,function(wcors_hc) {
    sapply(1:cutk, function(i){
      pairs = names(wcors_hc[wcors_hc == i])
      genes = sort(unique(unlist(strsplit(pairs,"_"))))
      return(genes)
    }, simplify = FALSE)
  })

  # perform pathway testing per cluster
  # want a list of dataframes, with column for cluster, pathway, p-value, numbergenes
  pathwayDF = sapply(names(sigWcorsList), function(branch) {
    dfBranch = sapply(1:length(clusterGenes[[branch]]),function(i){

      # print(i)
      # perfrom pathway test for this set of genes
      genes = clusterGenes[[branch]][[i]]

      hallmarks_pval = unlist(lapply(hallmarks, function(hallmark){
        if (!any(geneUniverse %in% hallmark)) return(1)
        if (!any(geneUniverse %in% genes)) return(1)
        fisher.test(table(geneUniverse %in% hallmark,geneUniverse %in% genes), alt = "g")$p.value
      }))

      df = data.frame(cluster = i,
                      pathway = names(hallmarks_pval),
                      pvalue = hallmarks_pval)
      return(df)
    }, simplify = FALSE)

    return(do.call(rbind,dfBranch))

  }, simplify = FALSE)

  message("Performed pathway analysis for branches")

  # calculate average wcor per cluster
  wcorAverage = sapply(names(sigWcorsList),function(branch) {

    t(sapply(unique(sig_wcors_hc[[branch]]), function(i){

      colMeans(sigWcorsList[[branch]][sig_wcors_hc[[branch]] == i,, drop = FALSE])

    }, simplify = TRUE))

  }, simplify = FALSE)

  message("Built pathway and average weighted correlation data frame")

  pathwaysToDisplay = lapply(pathwayDF, function(pathwaydf) {
    unique(unlist(sapply(unique(pathwaydf$cluster), function(i){
      # print(i)
      subdf = subset(pathwaydf, cluster == i)
      subdfp = subdf[,"pvalue"]
      names(subdfp) <- subdf[,"pathway"]
      names(sort(subdfp)[1:topPathways])
    }, simplify = FALSE)))
  })

  message("Extracted top pathways to display")


  pathwayDFsub = sapply(names(sigWcorsList), function(branch) {
    subset(pathwayDF[[branch]], pathway %in% pathwaysToDisplay[[branch]])
  }, simplify = FALSE)

  if (cluster) {

    require(ComplexHeatmap)

    hList = sapply(names(sigWcorsList), function(branch) {
      mat = cast(pathwayDFsub[[branch]], cluster ~ pathway)[,-1]
      h = Heatmap(-log10(mat),
                  cluster_rows = TRUE,
                  row_names_gp = gpar(fontsize = 8),
                  col = c("white","red","red"),
                  # name = "-log(P-value)",
                  heatmap_legend_param = list(legend_direction = "horizontal"),
                  rect_gp = gpar(col = "grey", lty = 1, lwd = 1),
                  row_names_max_width = unit(200,"mm")
      )
      return(h)
    })
    message("Completed clustering of weighted correlation types")
  }

  message("Extracted pathway clusters")

  pathway_gList = sapply(names(pathwayDFsub), function(branch) {

    pdf = pathwayDFsub[[branch]]
    # wrap the labels

    if (cluster) {
      pdf$cluster = factor(pdf$cluster, levels = row_order(hList[[branch]])[[1]])
      pdf$pathway = factor(pdf$pathway, levels = unique(pdf$pathway)[column_order(hList[[branch]])])
    }

    pdf$newpathway = factor(stringr::str_wrap(gsub("_"," ",pdf$pathway), width = label_wrap_cut),
                            levels = stringr::str_wrap(gsub("_"," ",levels(pdf$pathway)), width = label_wrap_cut))


    g = ggplot(pdf,aes(x = cluster, y  = newpathway,
                       fill = -log10(pvalue), size = -log10(pvalue),
                       color = -log10(pvalue))) +
      geom_point() +
      theme_minimal() +
      scale_fill_gradient(low = "grey",high = "red") +
      scale_color_gradient(low = "grey",high = "red") +
      xlab("") +
      ylab("") +
      NULL

    return(g)
  }, simplify = FALSE)

  message("Built pathway graphs")

  # weighted correlation graph
  wcor_gList = sapply(names(wcorAverage), function(branch){

    groupmeansdf = melt(wcorAverage[[branch]])
    colnames(groupmeansdf) <- c("X1","X2","value")

    if (cluster) {
      groupmeansdf$X1 = factor(groupmeansdf$X1,levels = row_order(hList[[branch]])[[1]])
      h_split = cutree(as.hclust(row_dend(hList[[branch]])[[1]]),min(6, cutk))
      groupmeansdf$split = factor(h_split[as.character(groupmeansdf$X1)])
    } else {
      groupmeansdf$X1 = factor(groupmeansdf$X1)
      groupmeansdf$split = groupmeansdf$X1
    }


    g = ggplot(groupmeansdf,
               aes(x = X2, y = value, group = X1, col = split)) +
      geom_line(size = 2, show.legend = FALSE) +
      theme_minimal() +
      geom_hline(yintercept = 0) +
      xlab("Pseudotime") +
      ylab("Local correlation") +
      facet_grid(~X1) +
      NULL

    return(g)
  }, simplify = FALSE)

  message("Built average weighted correlation graphs")

  patch_gList = sapply(names(pathwayDFsub),function(branch){
    pathway_gList[[branch]] + ggtitle(branch) +
      wcor_gList[[branch]] + theme(axis.text.x = element_blank()) +
      plot_layout(ncol = 1, heights = c(2,1))
  }, simplify = FALSE)

  message("combined all graphs")

  final_g = patchwork::wrap_plots(patch_gList, nrow = 1)
  return(final_g)
}
