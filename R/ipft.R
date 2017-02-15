## Rcpp::sourceCpp(normalizePath(file.path(".", "src", "ipf.cpp")))

setClassUnion("numericOrNULL",members=c("numeric", "NULL"))
setClassUnion("matrixOrNULL",members=c("matrix", "NULL"))

rmdata <- setClass(
  "rmdata",

  slots = c(
    NP  = "numeric",
    NG  = "numeric",
    ADP = "numeric",
    ADG = "numeric",
    AAD = "numeric",
    W   = "numeric",
    H   = "numeric",
    ANW = "numeric",
    AIW = "numeric",
    MIW = "numeric",
    SDIW = "numeric",
    SDNW = "numeric",

    anppc = "numeric",
    angpc = "numeric",
    sdnppc = "numeric",
    sdngpc = "numeric",
    anwpc = "numeric",
    sdnwpc = "numeric",
    aiwpc = "numeric",
    sdiwpc = "numeric",
    area = "numeric"

  )
)

ipfModel <- setClass(
  "ipfModel",

  slots = c(
    neighbours = "matrix",
    weights = "matrix",
    distances = "matrix",
    k = "numeric",
    groups = "numeric"
  )
)

ipfEstimation <- setClass(
  "ipfEstimation",

  slots = c(
    location = "matrix",
    grouploc = "matrixOrNULL",
    errors = "numericOrNULL"
  )
)

GRIDSIZE <- 5

#' Distance function
#'
#' This function computes the distance from every observation
#' in the test set to every observation in the train test
#'
#' @param train a vector, matrix or data frame containing a set
#'     of training examples
#' @param test a vector, matrix or data frame containing a set
#'     of test examples
#' @param method The method to be used to calculate the distance.
#'     Implemented methods are: 'euclidean', 'manhattan', 'norm',
#'     'LGD' and 'PLGD'
#' @param subset columns to use to compute the distance.
#' @param norm parameter for the 'norm' method
#' @param sd parameter for 'LGD' and 'PLGD' methods
#' @param epsilon parameter for 'LGD' and 'PLGD' methods
#' @param alpha parameter for 'PLGD' method
#' @param threshold parameter for 'PLGD' method
#'
#' @return This function returns a matrix with dimensions:
#'     nrow(test) x nrow(train), containing the distances from
#'     test observations to train observations
#'
#' @examples
#'     ipfDist(ipftrain[,1:168], ipftest[,1:168])
#'
#'     ipfDist(ipftrain, ipftest, subset = seq(1,168))
#'
#'     ipfDist(ipftrain, ipftest, subset = c('LONGITUDE', 'LATITUDE'), method = 'manhattan')
#'
#' @useDynLib ipft
#' @importFrom Rcpp sourceCpp
#'
#' @export
ipfDist <- function(train, test, method = 'euclidean', subset = NULL, norm = 2,
                    sd = 10, epsilon = 1e-30, alpha = 20, threshold = 20) {
  if (is.null(subset)) {
    m_train <- as.matrix(train)
    m_test <- as.matrix(test)
  } else {
    m_train <- as.matrix(train[, subset])
    m_test <- as.matrix(test[, subset])
  }
  if (ncol(m_train) == 1) {m_train <- t(m_train)}
  if (ncol(m_test) == 1) {m_test <- t(m_test)}
  switch(method,
         'euclidean' = dist <- ipfEuclidean(m_train, m_test),
         'manhattan' = dist <- ipfManhattan(m_train, m_test),
         'norm' = dist <- ipfNormDistance(m_train, m_test, norm),
         'LGD' = dist <- ipfLGD(m_train, m_test, sd, epsilon),
         'PLGD' = dist <- ipfPLGD(m_train, m_test, sd, epsilon, alpha, threshold),
         stop("invalid distance method.")
         )
  return (dist)
}

#' Transform function
#'
#' This function transforms the RSSI data to positive
#' or exponential values
#'
#' @param data    a vector, matrix or data frame containing the RSSI vectors
#' @param trans   the transformations to perform
#' @param minRSSI the minimum value for RSSI to consider when transforming
#'                the RSSI to positive values.
#' @param maxRSSI the maximum value for RSSI to consider when transforming
#'                the RSSI to exponential values.
#' @param noRSSI  value used in the RSSI data to represent a not detected AP.
#' @param alpha   parameter for exponential transformation
#'
#' @return This function returns a vector, matrix or data frame containing
#'         the transformed data
#'
#' @examples
#'     trainRSSI <- ipftrain[,1:168]
#'     ipfTransform(trainRSSI, trans = 'positive')
#'
#'     trainRSSI <- ipftrain[,1:168]
#'     posTrainRSSI <- ipfTransform(trainRSSI, trans = 'positive')
#'     expTrainRSSI <- ipfTransform(posTrainRSSI, trans = 'exponential', maxRSSI = 104)
#'
#' @export
ipfTransform <- function(data, trans = 'positive', minRSSI = -104, maxRSSI = 0,
                         noRSSI = 0, alpha = 24) {
  if (is.element('positive', trans)) {
    if (minRSSI >= 0) {stop("minRSSI must be negative.")}
    pdata <- data
    pdata[data == noRSSI] <- 0
    pdata[data != noRSSI] <- data[data != noRSSI] - minRSSI
    data <- pdata
  } else if (is.element('exponential', trans)) {
    data[data != noRSSI] <- exp(data[data != noRSSI] / alpha) / exp(maxRSSI / alpha)
  }

  return (data)
}

#' Creates groups based on the specified parameters
#'
#' This function groups the data based on the specified variables
#' and assigns an id to each group
#'
#' @param data A data frame
#' @param ...  Variables to group by. All variables (columns) will
#'             be used if no parameter is provided.
#'
#' @return     A numeric vector with the ids of the groups, in the
#'             same order as they appear in the data provided.
#'
#' @examples
#'
#'     ipfGroup(mtcars, cyl)
#'
#'     ipfGroup(mtcars, gear, carb)
#'
#'     ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'
#' @importFrom dplyr group_indices group_indices_
#'
#' @export
ipfGroup <- function(data, ...) {
  if (missing(..0)) {
    return (group_indices_(data, .dots = colnames(data)))
  } else {
    return (group_indices(data, ...))
  }
}

#' Creates clusters using the specified method
#'
#' This function creates clusters using the the specified method
#' and assigns a cluster id to each cluster
#'
#' @param data   a data frame
#' @param method the method to use to clusterize the data. Implemented
#'               methods are:
#'               'k-means' for k-means algorithm. Requires parameter k.
#'               'AP' for affinity propagation algorithm.
#' @param k      parameter k
#' @param ...    additional parameters for apcluster and apclusterK
#'
#' @return       A list with:
#'                  clusters -> a numeric vector with the ids of the clusters
#'                  centers  -> a data frame with the centers of the clusters
#'
#' @examples
#'
#'     ipfCluster(head(ipftrain, 20)[, 169:170], k = 4)
#'
#'     ipfCluster(head(ipftrain[, grep('^wap', names(ipftrain))], 20), method = 'AP')$clusters
#'
#' @importFrom stats kmeans
#' @importFrom dplyr select
#' @importFrom cluster pam
#' @importFrom apcluster apclusterK apcluster negDistMat
#'
#' @export
ipfCluster <- function(data, method = 'k-means', k = NULL, ...) {
  result <- NULL

  if (method == 'k-means') {
    if (!is.null(k)) {
      kmn <- kmeans(data, k)
    } else {
      stop("You must privide a value for k.")
    }
    result <- list(clusters = as.numeric(kmn$cluster), centers = data.frame(kmn$centers))
  } else if (method == 'k-medoid') {
    print("TODO: IMPEMENT k_MEDOID")
    result <- NULL
  } else if (method == 'AP') {
    if (!missing(k)) {
      ap <- apclusterK(negDistMat(r = 2), data, K = k, verbose = FALSE, ...)
    } else {
      ap <- apcluster(negDistMat(r = 2), data, ...)
    }
    result <- list(clusters = as.numeric(as.factor(ap@idx)),
                   centers = data.frame(data[ap@exemplars,]))
  } else {
    stop("invalid method selected.")
  }
  return (result)
}

#' This function implements the k-nearest neighbours algorithm
#'
#' @param train     a data frame containing the RSSI vectors of the training set
#' @param test      a data frame containing the RSSI vectors of the test set
#' @param k         the k parameter for knn algorithm (number of nearest neighbours)
#' @param method    the method to compute the distance between the RSSI vectors:
#'                 'euclidean', 'manhattan', 'norm', 'LGD' or 'PLGD'
#' @param norm      parameter for the 'norm' method
#' @param sd        parameter for 'LGD' and 'PLGD' methods
#' @param epsilon   parameter for 'LGD' and 'PLGD' methods
#' @param alpha     parameter for 'PLGD' method
#' @param threshold parameter for 'PLGD' method
#' @param FUN       an alternative function provided to compute the distance.
#'                  This function must return a matrix of dimensions:
#'                  nrow(test) x nrow(train), containing the distances from
#'                  test observations to train observations
#' @param ...       additional parameters for provided function FUN
#'
#' @return          An S4 class object of type ipfModel, with the following slots:
#'                  neighbours -> a matrix with k columns and nrow(test) rows, with the
#'                                k nearest neighbours for each test observation
#'                  weights ->    a matrix with k columns and nrow(test) rows, with the
#'                                weight for each neighbour
#'                  distances ->  a matrix with k columns and nrow(test) rows, with the
#'                                distances between test and each neighbour
#'                  k ->          k parameter
#'                  groups ->     the group index for each training observation
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168])
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168], k = 9, method = 'manhattan')
#'
#' @importFrom methods new
#'
#' @export
ipfKnn <- function(train, test, k = 3, method = 'euclidean', norm = 2, sd = 5,
                   epsilon = 1e-3, alpha = 1, threshold = 20, FUN = NULL, ...) {
  nr <- nrow(train)
  nt <- nrow(test)
  if (nr < k) stop("k can not be greater than the number of rows in training set")

  if (is.atomic(test)) {
    test <- matrix(test, nrow = 1)
  } else {
    test <- as.matrix(test)
  }
  train <- as.matrix(train)
  if (is.function(FUN)) {
    dm <- FUN(...)
  } else {
    dm <- ipfDist(train = train, test = test, method = method, norm = norm, sd = sd,
                  epsilon = epsilon, alpha = alpha, threshold = threshold)
  }
  ne <- matrix(0, nrow = nt, ncol = k)
  ws <- matrix(0, nrow = nt, ncol = k)
  for (i in 1:nt) {
    ne[i,] <- order(dm[i, ])[1:k]
    w <- (1 / (1 + dm[i, ne[i,]]))
    ws[i,] <- w / sum(w)
  }
  return(ipfModel(neighbours = ne, weights = ws, distances = dm, k = k, groups = seq(1, nr)))
}

#' This function implements a probabilistic algorithm
#'
#' @param train     a data frame containing the RSSI vectors of the training set
#' @param test      a data frame containing the RSSI vectors of the test set
#' @param groups    a numeric vector of length = nrow(train) containing the group index
#'                  for the training vectors
#' @param k         the k parameter for the algorithm (number of similar neighbours)
#' @param FUN       function to compute the similarity measurement. Default is 'sum'
#' @param delta     parameter delta
#' @param ...       additional parameters for provided function FUN
#'
#' @importFrom dplyr group_by summarise_each arrange "%>%" funs
#'
#' @return          An S4 class object of type ipfModel, with the following slots:
#'                  neighbours -> a matrix with k columns and nrow(test) rows, with the
#'                                k most similar training observation for each test observation
#'                  weights ->    a matrix with k columns and nrow(test) rows, with the weights
#'                  distances ->  a matrix with k columns and nrow(test) rows, with the distances
#'                  k ->          k parameter
#'                  groups ->     the group index for each training observation
#'
#' @examples
#'
#'     groups <- ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'     model <- ipfProb(ipftrain[, 1:168], ipftest[, 1:168], groups)
#'
#'     groups <- ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'     model <- ipfProb(ipftrain[, 1:168], ipftest[, 1:168], groups, k = 9, delta = 10)
#'
#' @importFrom methods new
#'
#' @export
ipfProb <- function(train, test, groups, k = 3, FUN = sum, delta = 1, ...) {
  if (ncol(train) != ncol(test)) {
    stop("Training and test datasets must have the same number of columns.")
  }
  if (nrow(train) != length(groups)) {
    stop("Incorrect dimmension of groups. Should be equal to nrow(train).")
  }

  GROUP <- NULL
  sd <- NULL

  train$GROUP <- groups
  tr_mn <- train %>% group_by(GROUP) %>% summarise_each(funs(mean)) %>% arrange(GROUP)
  tr_sd <- train %>% group_by(GROUP) %>% summarise_each(funs(sd)) %>% arrange(GROUP)

  tr_mn$GROUP <- NULL
  tr_sd$GROUP <- NULL

  nr <- nrow(tr_mn)
  nc <- ncol(tr_mn)
  if (nr < k) stop("k can not be greater than the number of rows in training set")
  sim <- matrix(0, nrow(test), nr)
  p <- matrix(0, nrow(test), nc)
  for (j in 1:nr) {
    m <- as.numeric(tr_mn[j, ])
    s <- as.numeric(tr_sd[j, ])

    for (i in 1:nc) {
      o1 <- test[, i] - delta
      o2 <- test[, i] + delta

      p1 <- pnorm(o1, mean = m[i], sd = s[i], lower.tail = FALSE)
      p2 <- pnorm(o2, mean = m[i], sd = s[i], lower.tail = FALSE)

      p[,i] <- (p1 - p2)
    }
    sim[, j] <- apply(p, 1, FUN, ...)
  }

  ne <- matrix(0, nrow = nrow(test), ncol = k)
  ws <- matrix(0, nrow = nrow(test), ncol = k)
  for (i in 1:nrow(test)) {
    ne[i,] <- order(sim[i, ], decreasing = TRUE)[1:k]
    w <- sim[i, ne[i,]]
    ws[i,] <- w / sum(w)
  }
  return(ipfModel(neighbours = ne, weights = ws, k = k, groups = groups))
}

#' This function estimates the location of the test observations
#'
#' @param ipfmodel an ipfModel
#' @param locdata  a matrix or a data frame containing the position of the training set observations
#' @param loctest  a matrix or a data frame containing the position of the test set observations
#'
#' @return         An S4 class object of type ipfEstimation, with the following slots:
#'                 location ->   a matrix with the predicted locations
#'                 grouploc ->   a matrix with the location data for each group
#'                 errors   ->   a numeric vector with the errors
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168])
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'
#'     groups <- ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'     model <- ipfProb(ipftrain[, 1:168], ipftest[, 1:168], groups, k = 9, delta = 10)
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'
#' @importFrom dplyr group_by summarise_each arrange "%>%" funs
#'
#' @importFrom methods new
#'
#' @export
ipfEstimate <- function(ipfmodel, locdata, loctest = NULL) {
  if (class(ipfmodel) != 'ipfModel') {stop("Wrong parameter type for wfmodel.")}
  mloc <- matrix(0, nrow(ipfmodel@neighbours), 2)
  grouploc <- as.matrix(locdata)
  GROUP <- NULL    # To avoid 'no visible binding for global variable' NOTE

  if (!is.null(ipfmodel@groups)) {
    locdata$GROUP <- ipfmodel@groups
    locdata <- locdata %>% group_by(GROUP) %>% summarise_each(funs(mean)) %>% arrange(GROUP)
    grouploc <- as.matrix(locdata)
    locdata$GROUP <- NULL
    locdata <- as.data.frame(locdata)
  }

  mloc[, 1] <- rowSums(locdata[ipfmodel@neighbours, 1] * ipfmodel@weights)
  mloc[, 2] <- rowSums(locdata[ipfmodel@neighbours, 2] * ipfmodel@weights)
  errors <- c()
  if (!is.null(loctest)) {
    errors <- sqrt((mloc[, 1] - loctest[, 1])^2 + (mloc[, 2] - loctest[, 2])^2)
  }
  return(ipfEstimation(location = mloc, grouploc = grouploc, errors = errors))
}

# AVERAGE NUMBER OF WAPS PER OBSERVATION
rmce_getANW <- function(traindata) {
  traindata[traindata != 0] = 1
  return (mean(rowSums(traindata)))
}

# STANDARD DEVIATION OF THE NUMBER OF WAPS PER OBSERVATION
rmce_getSDNW <- function(traindata) {
  traindata[traindata != 0] = 1
  return (sd(rowSums(traindata)))
}

# AVERAGE RSSI (WHEN > 0)
rmce_getAIW <- function(traindata) {
  return(mean(c(t(traindata[traindata != 0]))))
}

# MEDIAN RSSI (WHEN > 0)
rmce_getMIW <- function(traindata) {
  return(median(c(t(traindata[traindata != 0]))))
}

# STANDARD DEVIATION OF RSSI (WHEN > 0)
rmce_getSDIW <- function(traindata) {
  return(sd(c(t(traindata[traindata != 0]))))
}

# NUMBER OF OBSERVATIONS
rmce_getNP <- function(locdata) {
  return (nrow(locdata))
}

# NUMBER OF GROUPS
#' @importFrom dplyr group_indices_
rmce_getNG <- function(locdata) {
  # print(dim(locdata))
  # print(length(unique(group_indices_(locdata, .dots = names(locdata)))))
  return (length(unique(group_indices_(locdata, .dots = names(locdata)))))
}

# AVERAGE DENSITY OF OBSERVATIONS
rmce_getADP <- function(locdata) {
  return (mean(rowMeans(ipfDist(locdata, locdata))))
}

# AVERAGE DENSITY OF GROUPS
#' @importFrom dplyr group_indices_ group_by summarise_each arrange "%>%" funs
#'
rmce_getADG <- function(locdata) {
  GROUP <- NULL
  gr <- group_indices_(locdata, .dots = names(locdata))
  locdata$GROUP <- gr
  locdata <- locdata %>% group_by(GROUP) %>% summarise_each(funs(mean)) %>% arrange(GROUP)
  locdata$GROUP <- NULL
  locdata <- as.data.frame(locdata)
  return (rmce_getADP(locdata))
}

# AVERAGE DISTANCE BETWEEN TWO RANDOM POINTS INSIDE A RECTANGLE
rmce_getAAD <- function(locdata) {
  return (rectAD(rmce_getW(locdata), rmce_getH(locdata)))
}

# WIDTH
rmce_getW <- function(locdata) {
  minW <- min(locdata[, 1])
  maxW <- max(locdata[, 1])
  return (abs(maxW - minW))
}

# HEIGHT
rmce_getH <- function(locdata) {
  minH <- min(locdata[, 2])
  maxH <- max(locdata[, 2])
  return (abs(maxH - minH))
}

# TRAINING RADIO MAP DATA PER CELL
rm_getCD <- function(traindata, locdata, gridSize = 5) {
  minW <- min(locdata[, 1])
  maxW <- max(locdata[, 1])
  minH <- min(locdata[, 2])
  maxH <- max(locdata[, 2])

  columns <- ceiling(abs(maxW - minW) / gridSize)
  rows <- ceiling(abs(maxH - minH) / gridSize)

  nps <- c()
  ngs <- c()
  anw <- c()
  aiw <- c()
  miw <- c()
  sdw <- c()

  ng <- 0

  for (c in 1:columns) {
    for (r in 1:rows) {
      lxl <- minW + (c - 1) * gridSize
      lxr <- lxl + gridSize
      lyb <- minH + (r - 1) * gridSize
      lyt <- lyb + gridSize

      s <- locdata[, 1] >= lxl & locdata[, 1] < lxr & locdata[, 2] >= lyb & locdata[, 2] < lyt
      if (sum(s) == 0) {
        next
      }
      ng <- ng + 1
      p <- locdata[s, ]
      t <- traindata[s, ]

      nps <- c(nps, nrow(p))           # Number of observations in the grid
      ngs <- c(ngs, rmce_getNG(p))     # Number of groups in the grid
      anw <- c(anw, rmce_getANW(t))    # Average number of WAPs heard per observation in the grid
      aiw <- c(aiw, rmce_getAIW(t))    # Average RSSI in the grid
      miw <- c(miw, rmce_getMIW(t))    # Median RSSI in the grid
      sdw <- c(sdw, rmce_getSDIW(t))   # Stardard deviation of RSSI in the grid

    }
  }

  return (list(
    anppc = mean(nps), angpc = mean(ngs), sdnppc = sd(nps), sdngpc = sd(ngs),
    anwpc = mean(anw), sdnwpc = sd(anw), aiwpc = mean(aiw), sdiwpc = sd(aiw),
    area = ng * gridSize * gridSize
  ))
}

ipfrmdata <- function(traindata, locdata, gridSize = GRIDSIZE) {
  ANW  <- rmce_getANW(traindata)    # AVERAGE NUMBER OF WAPS PER OBSERVATION
  AIW  <- rmce_getAIW(traindata)    # AVERAGE RSSI
  MIW  <- rmce_getMIW(traindata)    # MEDIAN RSSI
  SDIW  <- rmce_getSDIW(traindata)  # STANDARD DEVIATION OF RSSI
  SDNW <- rmce_getSDNW(traindata)   # STANDARD DEVIATION OF RSSI

  NP <- rmce_getNP(locdata)     # TOTAL NUMBER OF OBSERVATIONS
  NG <- rmce_getNG(locdata)     # TOTAL NUMBER OF GROUPS
  ADP <- rmce_getADP(locdata)   # AVERAGE DISTANCE BETWEEN OBSERVATIONS
  ADG <- rmce_getADG(locdata)   # AVERAGE DISTANCE BETWEEN GROUPS
  AAD <- rmce_getAAD(locdata)   # AVERAGE DISTANCE BETWEEN TWO RANDOM POINTS (RECTANGLE WxH)
  W <- rmce_getW(locdata)       # TOTAL AREA WIDTH
  H <- rmce_getH(locdata)       # TOTAL AREA HEIGHT

  rmcd <- rm_getCD(traindata, locdata)

  RmData <- NULL

  wfrmdata <- RmData(
    ANW = ANW, AIW = AIW, MIW = MIW, SDIW = SDIW, SDNW = SDNW,
    NP = NP, NG = NG, ADP = ADP, ADG = ADG, AAD = AAD,
    W = W, H = H, anppc = rmcd$anppc, angpc = rmcd$angpc,
    sdnppc = rmcd$sdnppc, sdngpc = rmcd$sdngpc, anwpc = rmcd$anwpc,
    sdnwpc = rmcd$sdnwpc, aiwpc = rmcd$aiwpc, sdiwpc = rmcd$sdiwpc,
    area = rmcd$area
  )
  return (wfrmdata)
}

rmce <- function(rmdata) {
  NP <- rmdata@NP   # NUMBER OF POINTS
  NG <- rmdata@NG   # NUMBER OF GROUPS
  ADP <- rmdata@ADP # AVERAGE DISTANCE BETWEEN TRAINING POINTS
  ADG <- rmdata@ADG # AVERAGE DISTANCE BETWEEN TRAINING GROUPS
  AAD <- rmdata@AAD # ANALYTIC AVERAGE DISTANCE (RECTANGLE)
  W <- rmdata@W     # WIDTH
  H <- rmdata@H     # HEIGHT
  ANW <- rmdata@ANW # AVERAGE NUMBER OF WAPS PER OBSERVATION
  AIW <- rmdata@AIW # AVERAGE RSSI (WHEN > 0)
  MIW <- rmdata@MIW # MEDIAN RSSI (WHEN > 0)
  SDW <- rmdata@SDW # STANDARD DEVIATION OF RSSI (WHEN > 0)

  A <- W * H        # AREA

  DISP <- NG / NP   # DISPERSION OF POINTS

  # ce <- A/(N*K)
  # ce <- sqrt(G * N * A / MP)
  # ce <- MP / DISP / 100000
  # ce <- (ADP * AAD / DISP) / sqrt(A) / 50

  ce <- (AIW / 25) * (ADP * ADG * AAD) / (W * H)

  # print(paste0('ANW: ', round(ANW, 2), ', AIW: ', round(AIW, 2), ', MIW: ',
  #              round(MIW, 2), ', SDW: ', round(SDW, 2), ', NP: ', NP, ', NG: ',
  #              NG, ', ADP: ', round(ADP, 2), ', ADG: ', round(ADG, 2), ', AAD: ',
  #              round(AAD, 2), ', W: ', round(W, 2), ', H: ', round(H, 2), ', ce: ',
  #              round(ce, 2)))
  return (ce)
}

rectAD <- function(L1, L2) {
  d <- sqrt(L1^2 + L2^2)
  t1 <- (L1^3) / (L2^2) + (L2^3) / (L1^2)
  t2 <- d * (3 - (L1^2) / (L2^2) - (L2^2) / (L1^2))
  t3 <- (5 / 2) * (((L2^2) / L1) * log((L1 + d) / L2) + ((L1^2) / L2) * log((L2 + d) / L1))
  return ((t1 + t2 + t3) / 15)
}

#' Plots the probability density function of the estimated error
#'
#' @param estimation  an ipfEstimation
#' @param xlab        x-axis label
#' @param ylab        y-axis label
#' @param title       plot title
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168])
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'     ipfPlotPdf(estimation)
#'
#'     groups <- ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'     model <- ipfProb(ipftrain[, 1:168], ipftest[, 1:168], groups, k = 9, delta = 10)
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'     ipfPlotPdf(estimation, title = 'Probability density function')
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline scale_color_manual
#' @importFrom ggplot2 element_rect scale_alpha_manual scale_x_continuous labs theme
#' @importFrom methods is
#'
#' @export
ipfPlotPdf <- function(estimation, xlab = 'error', ylab = 'density', title = '') {
  if (isS4(estimation) && is(estimation, "ipfEstimation")) {
    if (is.null(estimation@errors)) {
      stop("The estimation has no errors data.")
    }
  } else {
    stop("ipfPlotPdf: Wrong parameter type. ipfEstimation expected.")
  }

  # Avoid no visible binding for global variable
  ..density.. <- NULL
  median <- NULL
  errors <- data.frame(ERRORS = estimation@errors)

  tickroundup <- c(1, 2, 5, 10)
  xmax <- max(errors)
  xmark <- xmax / 20
  xtick <- 10^floor(log10(xmark)) *
    tickroundup[[which(xmark <= 10^floor(log10(xmark)) * tickroundup)[[1]]]]

  ggplot(data = errors, aes(errors$ERRORS)) +
    geom_histogram(aes(y = ..density..), binwidth = max(errors$ERRORS) / 50, alpha = .5) +
    geom_density(alpha = .2, fill = "#666666", colour = 'grey30') +
    geom_vline(aes(xintercept = median(errors$ERRORS), color = "Median"),
               alpha = 0.5, linetype = "longdash", size = 1, show.legend = TRUE) +
    geom_vline(aes(xintercept = mean(errors$ERRORS), color = "Mean"), alpha = 0.5,
               linetype = "longdash", size = 1, show.legend = TRUE) +
    scale_color_manual(name = ' ', values = c(Median = "royalblue1", Mean = "firebrick1")) +
    scale_alpha_manual(values = c(0.5, 0.5)) +
    scale_x_continuous(breaks = seq(0, xmax, xtick)) +
    labs(x = "error",y = "density")  +
    labs(title = title) +
    theme(legend.position = c(0.9, 0.85),
          legend.background = element_rect(color = "transparent",
                                           fill = "transparent",
                                           size = 0,
                                           linetype = "blank"))
}

#' Plots the cumulative distribution function of the estimated error
#'
#' @param estimation  an ipfEstimation
#' @param xlab        x-axis label
#' @param ylab        y-axis label
#' @param title       plot title
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168])
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'     ipfPlotEcdf(estimation)
#'
#'     groups <- ipfGroup(ipftrain, LONGITUDE, LATITUDE)
#'     model <- ipfProb(ipftrain[, 1:168], ipftest[, 1:168], groups, k = 9, delta = 10)
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'     ipfPlotEcdf(estimation, title = 'Error cumulative distribution function')
#'
#' @importFrom ggplot2 ggplot aes stat_ecdf geom_vline scale_color_manual
#' @importFrom ggplot2 element_rect scale_alpha_manual scale_x_continuous labs theme
#' @importFrom methods is
#' @importFrom stats median pnorm sd
#'
#' @export
ipfPlotEcdf <- function(estimation, xlab = 'error',
                        ylab = 'cumulative density of error', title = '') {
  if (isS4(estimation) && is(estimation, "ipfEstimation")) {
    if (is.null(estimation@errors)) {
      stop("The estimation has no errors data.")
    }
  } else {
    stop("ipfPlotEcdf: Wrong parameter type. ipfEstimation expected.")
  }
  errors <- data.frame(ERRORS = estimation@errors)

  tickroundup <- c(1, 2, 5, 10)
  xmax <- max(errors)
  xmark <- xmax / 20
  xtick <- 10^floor(log10(xmark)) *
    tickroundup[[which(xmark <= 10^floor(log10(xmark)) * tickroundup)[[1]]]]

  ggplot(errors, aes(errors$ERRORS)) +
    stat_ecdf(geom = "point", alpha = .25) +
    geom_vline(aes(xintercept = median(errors$ERRORS), color="Median"), alpha = 0.5,
               linetype = "longdash", size = 1, show.legend = TRUE) +
    geom_vline(aes(xintercept = mean(errors$ERRORS), color="Mean"), alpha = 0.5,
               linetype = "longdash", size = 1, show.legend = TRUE) +
    scale_color_manual(name = ' ', values = c(Median = "royalblue1", Mean = "firebrick1")) +
    scale_alpha_manual(values = c(0.5, 0.5)) +
    scale_x_continuous(breaks = seq(0, xmax, xtick)) +
    labs(y = "cumulative density of error", x = "error") +
    labs(title = title) +
    theme(legend.position = c(0.9, 0.85),
          legend.background = element_rect(color = "transparent",
                                           fill = "transparent",
                                           size = 0,
                                           linetype = "blank"))
}

#' Plots the spatial location of the observations
#'
#' @param locdata     a data frame or matrix with the positions
#' @param plabel      if TRUE, adds labels to groups / observations
#' @param reverseAxis swaps axis
#' @param xlab        x-axis label
#' @param ylab        y-axis label
#' @param title       plot title
#' @param pgrid       plot grid (boolean)
#'
#' @examples
#'
#'     ipfPlotLoc(ipftrain[, 169:170])
#'
#'     ipfPlotLoc(ipftrain[, 169:170], plabel = TRUE, reverseAxis = TRUE,
#'                title = 'Position of training set observations')
#'
#' @importFrom dplyr group_by summarise_each arrange "%>%" funs
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_rect
#'
#' @export
ipfPlotLoc <- function(locdata, plabel = FALSE, reverseAxis = FALSE, xlab = NULL, ylab = NULL,
                       title = '', pgrid = FALSE) {

  if (reverseAxis) {
    locdata <- locdata[,c(2, 1)]
  }

  p <- ggplot(locdata)
  # To avoid no visible binding for global variable
  GROUP <- NULL
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL

  if (pgrid) {
    minX <- min(locdata[, 1])
    maxX <- max(locdata[, 1])
    minY <- min(locdata[, 2])
    maxY <- max(locdata[, 2])

    columns <- ceiling(abs(maxX - minX) / GRIDSIZE)
    rows <- ceiling(abs(maxY - minY) / GRIDSIZE)
    rect.data <- NULL

    for (c in 1:columns) {
      for (r in 1:rows) {
        lxl <- minX + (c - 1) * GRIDSIZE
        lxr <- lxl + GRIDSIZE
        lyb <- minY + (r - 1) * GRIDSIZE
        lyt <- lyb + GRIDSIZE

        s <- locdata[, 1] >= lxl & locdata[, 1] < lxr & locdata[, 2] >= lyb & locdata[, 2] < lyt

        if (sum(s) != 0) {
          rdata <- data.frame(xmin = lxl, xmax = lxr, ymin = lyb, ymax = lyt)
          rect.data <- rbind(rect.data, rdata)
        }
      }
    }
    suppressWarnings(
      p <- p + geom_rect(
        data = rect.data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, x = xmin, y = ymin),
        colour = 'grey60', alpha = 0.2, size = 0.3
      ) + geom_rect(
        aes(xmin = min(rect.data$xmin), xmax = max(rect.data$xmax),
            ymin = min(rect.data$ymin), ymax = max(rect.data$ymax),
            x = min(rect.data$xmin), y = min(rect.data$ymin)),
        colour = 'grey60', alpha = 0, size = 0.3, fill = 'transparent'
      )
    )


  }

  p <- p + geom_point(aes(x = locdata[, 1], y = locdata[, 2]), alpha = .25)

  if (plabel) {
    groups <- ipfGroup(locdata)
    label.data <- locdata
    label.data$GROUP = groups
    label.data <- label.data %>% group_by(GROUP) %>% summarise_each(funs(mean)) %>% arrange(GROUP)
    p <- p + geom_text(
      data = label.data,
      aes(x = label.data[, 2], y = label.data[, 3], label = label.data$GROUP),
      hjust = 0, vjust = 0, color = 'grey20', size = 3, nudge_x = 0.05, nudge_y = 0.05)
  }

  if (is.null(xlab)) {
    xlab <- colnames(locdata)[1]
  }

  if (is.null(ylab)) {
    ylab <- colnames(locdata)[2]
  }
  p <- p + labs(y = ylab, x = xlab) + labs(title = title)

  p

}

#' Plots the estimated locations
#'
#' @param model           an ipfModel
#' @param estimation      an ipfEstimation
#' @param testloc         location of test observations
#' @param observations    a numeric vector with the indices of estimations to plot
#' @param reverseAxis     swaps axis
#' @param showNeighbours  plot the k selected neighbours
#' @param showLabels      shows labels
#' @param xlab            x-axis label
#' @param ylab            y-axis label
#' @param title           plot title
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftest[, 1:168])
#'     estimation <- ipfEstimate(model, ipftrain[, 169:170], ipftest[, 169:170])
#'     ipfPlotEst(model, estimation, ipftest[, 169:170], observations = seq(7,10),
#'                showNeighbours = TRUE, reverseAxis = TRUE)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_curve arrow unit
#'
#' @export
ipfPlotEst <- function(model, estimation, testloc = NULL, observations = c(1),
                       reverseAxis = FALSE, showNeighbours = FALSE, showLabels = FALSE,
                       xlab = NULL, ylab = NULL, title = '') {

  ePoints <- as.data.frame(matrix(estimation@location[observations, ], length(observations), 2))
  if (reverseAxis) {
    if (!is.null(testloc)) {
      testloc <- testloc[,c(2, 1)]
    }
    ePoints <- ePoints[,c(2, 1)]
  }

  sPoints <- NULL

  p <- ggplot() + geom_point(data = ePoints, aes(x = ePoints[, 1],
                                                 y = ePoints[, 2]), size = 3, alpha = 0.5)

  if (showNeighbours) {
    nPoints <- NULL
    for (i in 1:length(observations)) {
      np <- as.data.frame(matrix(estimation@grouploc[model@neighbours[observations[i],], 2:3],
                                 model@k, 2))
      if (reverseAxis) {
        np <- np[,c(2, 1)]
      }
      nPoints <- rbind(nPoints, np)
      for (j in 1:model@k) {
        sPoints <- rbind(sPoints, cbind(ePoints[i, ], nPoints[j + (i - 1) * model@k, ]))
      }
    }
    p <- p + geom_point(data = nPoints, aes(x = nPoints[, 1], y = nPoints[, 2]),
                        size = 1.5, alpha = 0.75, col = 'orange') +
             geom_segment(data = sPoints, aes(x = sPoints[,1], xend = sPoints[, 3],
                                              y = sPoints[, 2], yend = sPoints[, 4]),
               alpha = 0.25, col = 'blue', size = 2)
  }

  if (!is.null(testloc)) {
    cPoints <- NULL
    for (i in 1:length(observations)) {
      cPoints <- rbind(cPoints, cbind(ePoints[i, ], testloc[observations[i], ]))
    }
    cPoints[cPoints[,1] == cPoints[,3] & cPoints[,2] == cPoints[,4], ] <- NA
    p <- p + geom_point(data = cPoints,
                       aes(x = cPoints[, 3], y = cPoints[, 4]),
                       size = 4, alpha = 0.5, col = 'green', na.rm = TRUE)
    p <- p + geom_curve(data = cPoints,
                        aes(x = cPoints[,3], xend = cPoints[, 1],
                            y = cPoints[, 4], yend = cPoints[, 2]),
                        alpha = 0.4, col = 'red', size = 1,
                        arrow = arrow(length = unit(0.03, "npc")),
                        curvature = 0.15, na.rm = TRUE)
  }

  if (showLabels) {
    ePoints$LABELS <- observations
    p <- p + geom_text(
      data = ePoints,
      aes(x = ePoints[, 1], y = ePoints[, 2], label = ePoints$LABELS),
      hjust = 0, vjust = 0, color = 'black', size = 4)
  }

  if (is.null(xlab)) {
    if (!is.null(testloc)) {
      xlab <- colnames(testloc)[1]
    } else {
      xlab <- ''
    }
  }

  if (is.null(ylab)) {
    if (!is.null(testloc)) {
     ylab <- colnames(testloc)[2]
    } else {
      ylab <- ''
    }
  }
  p <- p + labs(y = ylab, x = xlab) + labs(title = title)

  p
}
