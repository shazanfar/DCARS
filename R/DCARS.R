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

#' the DCARS function
#'
#' @title DCARS
#' @param dat a genes x samples gene expression rank matrix, should be already converted to ranks with first column lowest survival and last column highest survival
#' @param xname name of row of dat to test together with yname
#' @param yname name of row of dat to test together with xname
#' @param W weight matrix for weighted correlations,
#' @param rangeMin minimum range of weighted correlation vector to include for permutation testing
#' @param wcormin minimum absolute value weighted correlation vector to include for permutation testing
#' @param statmin minimum value DCARS test statistic to include for permutation testing
#' @param plot if TRUE plot observed weighted correlatin vector
#' @param niter number of iterations for permutation testing
#' @param extractTestStatisticOnly if TRUE, extract only the DCARS test statistic without permutation testing
#' @param extractWcorSequence if TRUE, extract only the weighted correlation vector without permutation testing
#' @param extractPermutationTestStatistics if TRUE, extract only the DCARS test statistic without permutation testing
#' @param verbose if TRUE, print updates
#' @param ... additional arguments passing on to weightMatrix()

#' @return either single value (p-value or test statistic), vector (local weighted correlation), or list (combination of above) depending on the input parameters

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


DCARS = function(dat,xname,yname,W=NULL,rangeMin = 0,wcormin = 0,statmin = 0,extractTestStatisticOnly = FALSE,extractWcorSequenceOnly = FALSE,plot = FALSE,niter = 100,extractPermutationTestStatistics = FALSE,verbose = FALSE,...) {

  # dat: is a genes x samples gene expression rank matrix, should be already converted to ranks with first column lowest survival and last column highest survival
  # xname: name of row of dat to test together with yname
  # yname: name of row of dat to test together with xname
  # W: weight matrix for weighted correlations,
  # rangeMin: minimum range of weighted correlation vector to include for permutation testing
  # wcormin: minimum absolute value weighted correlation vector to include for permutation testing
  # statmin: minimum value DCARS test statistic to include for permutation testing
  # plot: if TRUE plot observed weighted correlatin vector
  # niter: number of iterations for permutation testing
  # extractTestStatisticOnly: if TRUE, extract only the DCARS test statistic without permutation testing
  # extractWcorSequence: if TRUE, extract only the weighted correlation vector without permutation testing
  # extractPermutationTestStatistics if TRUE extracts the set of permuted DCARS statistics
  # verbose: if TRUE, print updates
  # ...: additional arguments passing on to weightMatrix()

  if (any(!c(xname,yname) %in% rownames(dat))) {
    message("xname and/or yname are not found in rownames(dat), returning NA")
    return(NA)
  }

  if (verbose) {
    print(paste(xname,yname))
  }

  if (is.null(W)) {
    message("weight matrix not specified, generating weight matrix now")
    W = weightMatrix(ncol(dat),...)
  }

  x = dat[xname,]
  y = dat[yname,]

  # weighted correlation vector
  wcor = sapply(1:length(x),function(i)weightedcor(x,y,W[i,]))
  if (any(!is.finite(wcor))) {
    i = which(!is.finite(wcor))
    wcor[i] <- mean(c(wcor[c(i-1,i+1)]))
  }

  if (plot) {
    plot(wcor,type="l",ylim=c(-1,1),lwd=3, col = "blue");abline(h=0,lty=2)
  }

  stat = sd(wcor)

  if (extractWcorSequenceOnly&extractTestStatisticOnly) {
    message("both extractWcorSequenceOnly and extractTestStatisticOnly have been specified as TRUE, returning both as a list")
    return(list(wcorSequence = wcor,
                TestStatistic = stat))
  }
  if (extractWcorSequenceOnly) return(wcor)
  if (extractTestStatisticOnly) return(stat)
  if (stat < statmin) return(NA)
  if (diff(range(wcor)) < rangeMin) return(NA)
  if (max(abs(wcor)) < wcormin) return(NA)

  # perform permutation test
  if (verbose) {
    message("performing permutation test")
  }
  sds = replicate(niter,{
    o = sample(1:length(x))
    xo = x[o]
    yo = y[o]
    wcor = sapply(1:length(x),function(i)weightedcor(xo,yo,W[i,]))
    if (any(!is.finite(wcor))) {
      i = which(!is.finite(wcor))
      wcor[i] <- mean(c(wcor[c(i-1,i+1)]))
    }
    sd(wcor)
  })

  if (extractPermutationTestStatistics) {
    if (verbose) {
      print("extract permuted test statistics")
    }
    # return(list(PermutedTestStatistics=sds))
    return(list(PermutedTestStatistics = sds))
  }

  return(mean(sds>=stat))
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
      edgeNames = apply(edgelist,1,function(x)paste0(sort(x[1:2]),collapse="_"))
    } else {
    edgeNames = rownames(edgelist)
    }
  }

  DCARSresult = apply(edgelist[,1:2],1,function(x)DCARS(dat=dat,xname = x[1],yname = x[2],...))
  if (is.matrix(DCARSresult)) {
    colnames(DCARSresult) <- edgeNames
  } else {
  names(DCARSresult) <- edgeNames
  }

  return(DCARSresult)
}

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


#' calculates weighted correlation between x and y
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

  # calculates weighted correlation between x and y
  # x and y are data vectors
  # w is a weight vector

  nw = sum(w)
  wssx = nw*sum(w*(x^2)) - sum(w*x)^2
  wssy = nw*sum(w*(y^2)) - sum(w*y)^2
  wssxy = nw*sum(w*x*y) - sum(w*x)*sum(w*y)
  wcor = wssxy/sqrt(wssx*wssy)
  return(wcor)
}

#' the stratifiedSample function
#'
#' @title stratifiedSample
#' @param stats a vector of statistics
#' @return \code{vector} a vector of indices for the statistics that should be used for permutation testing

#' @examples
#'
#'stats = rnorm(5000)
#'sampleindices = stratifiedSample(stats)
#'length(sampleindices)
#'
#' @export
#'
stratifiedSample = function(stats) {
  # takes in a vector of statistics, and returns a vector of indices for the statistics that should be used for permutation testing
  ranges = range(stats)
  fac = cut(stats,seq(from = 0, to = ranges[2], length.out = 101))
  nsamps = 3
  sampleindices = unlist(tapply(1:length(stats),fac,function(x){
    if (length(x)<nsamps) return(x)
    sample(x,nsamps,replace = FALSE)
  }))
  sampleindices = unique(sampleindices[!is.na(sampleindices)])
  # length(sampleindices)
  return(sampleindices)
}



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
