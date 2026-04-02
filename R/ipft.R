## Rcpp::sourceCpp(normalizePath(file.path(".", "src", "ipf.cpp")))

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
#'     dist <- ipfDistance(ipftrain[1:50, 1:20], ipftest[1:20, 1:20])
#'
#'     dist <- ipfDistance(ipftrain[1:50, ], ipftest[1:20, ], subset = c('X', 'Y'),
#'                        method = 'manhattan')
#'
#' @useDynLib ipft
#' @importFrom Rcpp sourceCpp
#'
#' @export
ipfDistance <- function(train, test, method = 'euclidean', subset = NULL, norm = 2,
                        sd = 10, epsilon = 1e-30, alpha = 20, threshold = 20) {
  train_has_dim <- !is.null(dim(train))
  test_has_dim <- !is.null(dim(test))

  if (is.null(subset)) {
    m_train <- as.matrix(train)
    m_test  <- as.matrix(test)
  } else {
    if (!train_has_dim || !test_has_dim) {
      stop("subset can only be used when train and test are matrices or data frames.")
    }
    m_train <- as.matrix(train[, subset, drop = FALSE])
    m_test  <- as.matrix(test[, subset, drop = FALSE])
  }

  # Vectors are treated as a single observation (one row).
  # One-column matrices/data.frames are treated as multiple observations (one column feature).
  if (!train_has_dim && ncol(m_train) == 1L) m_train <- t(m_train)
  if (!test_has_dim && ncol(m_test) == 1L) m_test  <- t(m_test)

  if (ncol(m_train) != ncol(m_test)) {
    stop("train and test must have the same number of columns.")
  }

  dist <- switch(method,
                 euclidean = ipfEuclidean(m_train, m_test),
                 manhattan = ipfManhattan(m_train, m_test),
                 norm      = ipfNormDistance(m_train, m_test, norm),
                 LGD       = ipfLGD(m_train, m_test, sd, epsilon),
                 PLGD      = ipfPLGD(m_train, m_test, sd, epsilon, alpha, threshold),
                 stop("invalid distance method.")
  )
  dist
}

#' Transform function
#'
#' Transforms the RSSI (Received Signal Strength Intensity) data to positive
#' or exponential values
#'
#' @param data      a vector, matrix or data frame containing the RSSI vectors
#' @param inRange   a vector containing the range of the RSSI value from the
#'                  initial data
#' @param outRange  the desired range for the output RSSI data.
#' @param inNoRSSI  value used in the RSSI data to represent a not detected AP.
#' @param outNoRSSI value desired in the RSSI output data to represent a not detected AP.
#' @param trans     the transformation to perform, 'scale' or 'exponential'
#' @param base      base for the 'exponential' transformation
#' @param alpha     alpha parameter for the 'exponential' transformation
#'
#' @return This function returns a vector, matrix or data frame containing
#'         the transformed data
#'
#' @examples
#'     trainRSSI <- ipftrain[,1:168]
#'     ipfTransform(trainRSSI, inRange = c(-100, 0), outRange = c(1, 100),
#'                  inNoRSSI = NA, outNoRSSI = 0)
#'
#' @export
ipfTransform <- function(data, outRange = c(0, 1), outNoRSSI = 0, inRange = NULL,
                         inNoRSSI = 0, trans = 'scale', base = exp(1), alpha = 24) {
  if (is.null(inRange)) {
    data_for_range <- as.matrix(data)
    if (!is.na(inNoRSSI)) {
      data_for_range[data_for_range == inNoRSSI] <- NA
    }
    mx <- suppressWarnings(max(data_for_range, na.rm = TRUE))
    mn <- suppressWarnings(min(data_for_range, na.rm = TRUE))
    if (!is.finite(mx) || !is.finite(mn)) {
      stop("ipfTransform: unable to infer inRange from data; provide inRange explicitly.")
    }
    inRange <- c(mn, mx)
  }

  if (!is.numeric(inRange) || length(inRange) != 2) {
    stop("ipfTransform: inRange must be a numeric vector of length 2.")
  }
  if (inRange[1] == inRange[2]) {
    stop("ipfTransform: inRange values must be different.")
  }
  if (!is.numeric(outRange) || length(outRange) != 2) {
    stop("ipfTransform: outRange must be a numeric vector of length 2.")
  }

  transform_fun <- NULL
  if (trans == 'scale') {
    transform_fun <- function(x) trans_lin(x, imin = inRange[1], imax = inRange[2],
                                           omin = outRange[1], omax = outRange[2],
                                           ins = inNoRSSI, ons = outNoRSSI)
  } else if (trans == 'exponential') {
    transform_fun <- function(x) trans_exp(x, imin = inRange[1], imax = inRange[2],
                                           omin = outRange[1], omax = outRange[2],
                                           ins = inNoRSSI, a = alpha, base = base)
  } else {
    stop("ipfTransform: invalid transformation.")
  }

  if (is.data.frame(data)) {
    tdata <- data
    tdata[] <- lapply(data, transform_fun)
    return(tdata)
  }

  transform_fun(data)
}

trans_lin <- function(x, imin, imax, omin, omax, ins,  ons) {
  if (is.na(ins)) {
    nos <- is.na(x)
  } else {
    nos <- x == ins
  }

  m <- (omin - omax) / (imin - imax)
  a <- omin - imin * m
  x <- a + x * m
  x[nos] <- ons
  x
}

trans_exp <- function(x, imin, imax, omin, omax, ins, a, base = exp(1)) {
  x <- trans_lin(x, imin, imax, omin = 0, omax = imax - imin, ins = ins, ons = 0)
  x <- base^(x/a) / base^(-imin/a)
  x
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
#'     group <- ipfGroup(mtcars, cyl)
#'
#'     group <- ipfGroup(mtcars, gear, carb)
#'
#'     group <- ipfGroup(ipftrain, X, Y)
#'
#' @importFrom dplyr group_indices
#'
#' @export
ipfGroup <- function(data, ...) {
  if (missing(..0)) {
    data <- as.data.frame(data)
    if (ncol(data) == 0) {
      return(rep(1L, nrow(data)))
    }
    grouped <- dplyr::group_by(data, dplyr::across(dplyr::all_of(colnames(data))))
    return(dplyr::group_indices(grouped))
  }
  grouped <- dplyr::group_by(as.data.frame(data), ...)
  dplyr::group_indices(grouped)
}

#' Creates clusters using the specified method
#'
#' Creates clusters using the the specified method
#' and assigns a cluster id to each cluster
#'
#' @param data   a data frame
#' @param method the method to use to clusterize the data. Implemented
#'               methods are:
#'               'k-means' for k-means algorithm. Requires parameter k.
#'               'grid' for clustering based on grid partition. Requires parameter grid.
#'               'AP' for affinity propagation algorithm.
#' @param k      parameter k
#' @param grid   a vector with the grid size for the 'grid' method
#' @param ...    additional parameters for k-means, apcluster and apclusterK
#'               for 'k-means' method additional parameters see:
#'               https://stat.ethz.ch/R-manual/R-devel/library/stats/html/kmeans.html
#'               for 'apcluster'AP' method additional parameters see:
#'               https://cran.r-project.org/web/packages/apcluster/index.html
#'
#' @return       A list with:
#'                  clusters -> a numeric vector with the ids of the clusters
#'                  centers  -> a data frame with the centers of the clusters
#'
#' @examples
#'
#'     clusters <- ipfCluster(head(ipftrain, 20)[, 169:170], k = 4)
#'
#'     clusters <- ipfCluster(head(ipftrain[, grep('^wap', names(ipftrain))], 20),
#'     method = 'AP')$clusters
#'
#' @importFrom stats kmeans
#' @importFrom dplyr select
#' @importFrom cluster pam
#' @importFrom apcluster apclusterK apcluster negDistMat
#'
#' @export
ipfCluster <- function(data, method = 'k-means', k = NULL, grid = NULL, ...) {
  result <- NULL

  if (method == 'k-means') {
    if (!is.null(k)) {
      kmn <- kmeans(x = data, centers = k, ...)
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
  } else if (method == 'grid') {
    if (is.null(grid)) stop("You must privide a value for grid.")
    if (ncol(data) != length(grid)) stop("The grid dimmensions must be equal to the data dimmensions.")
    nacs <- apply(data, 2, function(x) all(is.na(x)))
    data <- data[, !nacs]
    grid <- grid[!nacs]
    mins <- apply(data, 2, FUN = function(x) min(x, na.rm = TRUE))
    maxs <- apply(data, 2, FUN = function(x) max(x, na.rm = TRUE))
    dmns <- (maxs - mins) / grid

    dimg <- sapply(seq_along(grid), function(i) round(1 + (data[,i] - mins[i]) / grid[i]))

    dmgr <- length(grid)
    gmtx <- matrix(1, nrow(data), dmgr)

    for (i in 2:dmgr) {
      gmtx[,i] <- gmtx[,i - 1] * dmns[i]
    }

    groups <- as.numeric(as.factor(rowSums(dimg * gmtx)))

    centers <- matrix(0, length(unique(groups)), dmgr)

    for (i in unique(groups)) {
      centers[i, ] = colMeans(data[groups == i,])
    }

    centers <- data.frame(centers)
    colnames(centers) <- colnames(data)
    result <- list(clusters = groups, centers = centers)
  } else {
    stop("invalid method selected.")
  }
  return(result)
}

#' Estimates the positions of the emitter beacons
#'
#' @param fingerprints a data frame or a matrix with the RSSI fingerprints
#' @param positions    a data frame or a matrix with the positions of the fingerprints
#' @param method       method to use to estimate the position of the access points:
#'                     'centroid', 'wcentroid' or 'wip'
#' @param rssirange    a numeric vector with the range of the RSSI data
#' @param norssi       value used in dataRSSI when a beacon is not detected
#'
#' @examples
#'
#'     wapp <- ipfEstimateBeaconPositions(
#'       ipftrain[1:100, 1:168], ipftrain[1:100, 169:170],
#'       method = 'wcentroid'
#'     )
#'
#' @export
ipfEstimateBeaconPositions <- function(fingerprints, positions, method = 'wcentroid', rssirange = c(-100, 0),
                                       norssi = NA) {

  fp <- as.matrix(fingerprints)
  pos <- as.matrix(positions)

  if (nrow(fp) != nrow(pos)) {
    stop("fingerprints and positions must have the same number of rows.")
  }

  if (method == 'wcentroid') {
    tdata <- ipfTransform(fp, inRange = rssirange, outRange = c(0, 1),
                          inNoRSSI = norssi, outNoRSSI = 0)
  } else if (method == 'centroid') {
    tdata <- fp
    if (is.na(norssi)) {
      tdata[!is.na(tdata)] <- 1
      tdata[is.na(tdata)] <- 0
    } else {
      tdata[!is.na(tdata) & tdata != norssi] <- 1
      tdata[is.na(tdata) | tdata == norssi] <- 0
    }
  } else {
    stop("invalid method selected.")
  }

  # Normalize weights per beacon (column).
  col_sums <- colSums(tdata, na.rm = TRUE)
  tdata <- sweep(tdata, 2, col_sums, "/")
  tdata[, col_sums == 0] <- NA_real_

  aploc <- matrix(NA_real_, ncol(fp), ncol(pos))
  for (i in seq_len(ncol(fp))) {
    if (col_sums[i] == 0) next
    aploc[i, ] <- colSums(pos * tdata[, i], na.rm = TRUE)
  }
  aploc[is.nan(aploc)] <- NA
  aploc
}


#' Implements the k-nearest neighbors algorithm
#'
#' @param train_fgp a data frame containing the fingerprint vectors of the training set
#' @param train_pos a data frame containing the positions of the training set observations
#' @param k         the k parameter for knn algorithm (number of nearest neighbors)
#' @param method    the method to compute the distance between the RSSI vectors:
#'                 'euclidean', 'manhattan', 'norm', 'LGD' or 'PLGD'
#' @param weights   the algorithm to compute the weights: 'distance' or 'uniform'
#' @param norm      parameter for the 'norm' method
#' @param sd        parameter for 'LGD' and 'PLGD' methods
#' @param epsilon   parameter for 'LGD' and 'PLGD' methods
#' @param alpha     parameter for 'PLGD' method
#' @param threshold parameter for 'PLGD' method
#' @param FUN       an alternative function provided to compute the distance.
#'                  This function must return a matrix of dimensions:
#'                  nrow(test) x nrow(train), containing the distances from
#'                  test observations to train observations. The two first parameters
#'                  taken by the function must be train and test
#' @param ...       additional parameters for provided function FUN
#'
#' @return          An S3 object of class ipfModel, with the following properties:
#'                  params ->     a list with the parameters passed to the function
#'                  data   ->     a list with the fingerprints and locations
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[, 1:168], ipftrain[, 169:170], k = 9, method = 'manhattan')
#'
#' @importFrom methods new
#'
#' @export
ipfKnn <- function(train_fgp, train_pos, k = 3, method = 'euclidean', weights = 'distance',
                   norm = 2, sd = 5, epsilon = 1e-3, alpha = 1, threshold = 20, FUN = NULL, ...) {
  nr <- nrow(train_fgp)
  if (nr < k) stop("k can not be greater than the number of rows in training set")

  extra_params <- NULL
  if (!missing(..0)) {
    extra_params = list(...)
  }

  return(structure(list(
    params = list(name = 'knn',
                  k = k,
                  method = method,
                  norm = norm,
                  weights = weights,
                  sd = sd,
                  epsilon = epsilon,
                  alpha = alpha,
                  threshold = threshold,
                  FUN = FUN,
                  extra_params = extra_params),
    data   = list(fingerprints = as.data.frame(train_fgp),
                  positions = as.data.frame(train_pos))
  ), class = "ipfModel"))
}

knn <- function(tr_fgp, ts_fgp, k = 3, method = 'euclidean', weights = 'distance', norm = 2, sd = 5,
                epsilon = 1e-3, alpha = 1, threshold = 20, FUN = NULL, ...) {

  if (is.atomic(ts_fgp)) {
    ts_fgp <- matrix(ts_fgp, nrow = 1)
  } else {
    ts_fgp <- as.matrix(ts_fgp)
  }
  tr_fgp <- as.matrix(tr_fgp)
  if (is.function(FUN)) {
    if (missing(..0)) {
      dm <- FUN(tr_fgp, ts_fgp)
    } else {
      dm <- FUN(tr_fgp, ts_fgp, ...)
    }
  } else {
    dm <- ipfDistance(train = tr_fgp, test = ts_fgp, method = method, norm = norm, sd = sd,
                  epsilon = epsilon, alpha = alpha, threshold = threshold)
  }
  ne <- matrix(0, nrow = nrow(ts_fgp), ncol = k)
  ws <- matrix(0, nrow = nrow(ts_fgp), ncol = k)
  for (i in 1:nrow(ts_fgp)) {
    ne[i,] <- order(dm[i, ])[1:k]
    if (weights == 'distance') {
      w <- (1 / (1 + dm[i, ne[i,]]))
      ws[i,] <- w / sum(w)
    } else if (weights == 'uniform') {
      ws[i,] <- 1 / k
    }
  }
  return(list(neighbors = ne, weights = ws))
}

#' This function implements a probabilistic algorithm
#'
#' @param train_fgp  a data frame containing the fingerprint vectors of the training set
#' @param train_pos  a data frame containing the positions of the training set observations
#' @param group_cols a character vector with the names of the columns to be used as the criteria
#'                   to group the fingerprints. By default the groups will be created using all
#'                   the columns available in the train_pos data frame.
#' @param groups     a numeric vector of length = nrow(train) containing the group index
#'                   for the training vectors
#' @param k          the k parameter for the algorithm (number of similar neighbors)
#' @param FUN        function to compute the similarity measurement. Default is 'sum'
#' @param delta      parameter delta
#' @param ...        additional parameters for provided function FUN
#'
#' @importFrom dplyr group_by arrange "%>%"
#'
#' @return          An S3 object of class ipfModel, with the following properties:
#'                  params ->     a list with the parameters passed to the function
#'                  data   ->     a list with the fingerprints probabilistic parameters
#'                                (means and standard deviations) and its locations
#'
#' @examples
#'
#'     train <- ipftrain[1:200, ]
#'     groups <- ipfGroup(train, X, Y)
#'     model <- ipfProbabilistic(train[, 1:168], train[, 169:170], groups = groups)
#'
#'\dontrun{
#'     model <- ipfProbabilistic(ipftrain[, 1:168], ipftrain[, 169:170], k = 9, delta = 10)
#' }
#'
#' @importFrom methods new
#'
#' @export
ipfProbabilistic <- function(train_fgp, train_pos, group_cols = NULL, groups = NULL, k = 3,
                             FUN = sum, delta = 1, ...) {
  if (nrow(train_fgp) != nrow(train_pos)) {
    stop("train_fgp and train_pos must have the same number of rows")
  }
  if (!is.null(groups) && (nrow(train_fgp) != length(groups))) {
    stop("Incorrect dimmension of groups. Should be equal to the number of training fingerprints.")
  }

  train_fgp <- as.data.frame(train_fgp)
  train_pos <- as.data.frame(train_pos)

  GROUP <- NULL
  sd <- NULL
  meanNA <- function(x) mean(x, na.rm = TRUE)
  fPos <- function(x) if (!is.factor(x)) mean(x, na.rm = TRUE) else x[1]
  sdNA <- function(x) sd(x, na.rm = TRUE)

  if (is.null(groups)) {
    if (is.null(group_cols)) {
      groups <- ipfGroup(train_pos)
    } else {
      group_data <- cbind(train_fgp, train_pos)
      unknown_cols <- setdiff(group_cols, colnames(group_data))
      if (length(unknown_cols) > 0) {
        stop("Unknown columns in group_cols: ", paste(unknown_cols, collapse = ", "))
      }
      grouped <- dplyr::group_by(group_data, dplyr::across(dplyr::all_of(group_cols)))
      groups <- dplyr::group_indices(grouped)
    }
  }
  train_fgp$GROUP <- groups
  train_pos$GROUP <- groups
  tr_mns <- train_fgp %>%
    group_by(GROUP) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), meanNA), .groups = "drop") %>%
    arrange(GROUP)
  tr_sds <- train_fgp %>%
    group_by(GROUP) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sdNA), .groups = "drop") %>%
    arrange(GROUP)
  tr_pos <- train_pos %>%
    group_by(GROUP) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), fPos), .groups = "drop") %>%
    arrange(GROUP)

  tr_mns$GROUP <- NULL
  tr_sds$GROUP <- NULL
  tr_pos$GROUP <- NULL

  tr_mns[is.na(tr_mns)] <- NA  # NaN to NA

  extra_params <- NULL
  if (!missing(..0)) {
    extra_params = list(...)
  }

  return(structure(list(
    params = list(name = 'prob',
                  k = k,
                  FUN = FUN,
                  delta = delta,
                  extra_params = extra_params),
    data   = list(means = as.data.frame(tr_mns),
                  sds = as.data.frame(tr_sds),
                  positions = as.data.frame(tr_pos))
  ), class = "ipfModel"))

}

prob <- function(tr_mns, tr_sds, ts_fgp, k = 3, FUN = sum, delta = 1, ...) {

  nr <- nrow(tr_mns)
  nc <- ncol(tr_mns)
  if (nr < k) stop("k can not be greater than the number of rows in training set")
  sim <- matrix(0, nrow(ts_fgp), nr)
  p <- matrix(0, nrow(ts_fgp), nc)

  for (j in 1:nr) {
    m <- as.numeric(tr_mns[j, ])
    s <- as.numeric(tr_sds[j, ])

    for (i in 1:nc) {
      o1 <- ts_fgp[, i] - delta
      o2 <- ts_fgp[, i] + delta

      p1 <- pnorm(o1, mean = m[i], sd = s[i], lower.tail = FALSE)
      p2 <- pnorm(o2, mean = m[i], sd = s[i], lower.tail = FALSE)

      p[,i] <- (p1 - p2)
    }
    p[is.na(p)] <- 0
    sim[, j] <- apply(p, 1, FUN, ...)
  }
  ne <- matrix(0, nrow = nrow(ts_fgp), ncol = k)
  ws <- matrix(0, nrow = nrow(ts_fgp), ncol = k)
  for (i in 1:nrow(ts_fgp)) {
    ne[i,] <- order(sim[i, ], decreasing = TRUE)[1:k]
    w <- sim[i, ne[i,]]
    ws[i,] <- w / sum(w)
  }
  return(list(neighbors = ne, weights = ws))

}

#' Estimates the position of the observations from its fingerprints and the access point location
#' usins a logarithmic path loss model
#'
#' @param bpos         a matrix or a data frame containing the position of the beacons,
#'                     in the same order as they appear in fingerprints
#' @param rssirange    range of the RSSI data
#' @param norssi       value used to represent a not detected AP
#' @param alpha        path loss exponent
#' @param wapPow1      detected RSSI at one meter range
#'
#' @return             An S3 object of class ipfEstimation, with the following properties:
#'                     location  ->   a matrix with the predicted locations
#'                     errors    ->   a numeric vector with the errors, if loctest has been provided
#'                     neighbors ->   NULL
#'                     weights   ->   NULL
#'
#' @examples
#'
#'     ipfEst <- ipfProximity(ipftrain[1:10, 1:168], ipfpwap, ipftrain[1:10, 169:170], alpha = 4)
#'
#' @importFrom methods new
#'
#' @export
ipfProximity <- function(bpos, rssirange = c(-100, 0), norssi = NA, alpha = 5, wapPow1 = -30) {

  return(structure(list(
      params = list(name = 'prox',
                   rssirange = rssirange,
                   norssi = norssi,
                   alpha = alpha,
                   wapPow1 = wapPow1),
      data   = list(positions = as.data.frame(bpos))
  ), class = "ipfModel"))
}


prox <- function(fingerprints, fgps_pos, bpos, rssirange, norssi, alpha, wapPow1) {

  if (ncol(fingerprints) != nrow(bpos)) stop("Incorrect dimmensions.")
  notNAs <- !as.logical(rowSums(is.na(bpos[,])))
  fingerprints <- fingerprints[, notNAs]
  bpos <- bpos[notNAs, ]


  fingerprints[] <- lapply(fingerprints, ipfTransform, inRange = rssirange, outRange = c(0, 100),
                           inNoRSSI = norssi, outNoRSSI = 0)

  wP1 <- rep(trans_lin(x = wapPow1, imin = rssirange[1], imax = rssirange[2],
                       omin = 0, omax = 100, ins = norssi, ons = 0), ncol(fingerprints))

  mloc <- t(sapply(1:nrow(fingerprints), function(i) minim(fingerprint = fingerprints[i, ],
                                                           wapPos = bpos, wP1 = wP1,
                                                           alpha = alpha)))

  errors <- c()
  if (!is.null(fgps_pos)) {
    errors <- sqrt((mloc[, 1] - fgps_pos[, 1])^2 + (mloc[, 2] - fgps_pos[, 2])^2)
  }
  mloc <- as.data.frame(mloc)
  names(mloc) <- names(bpos)

  return(structure(list(
    location = mloc,
    errors = errors,
    confusion = NULL,
    neighbors = NULL,
    weights = NULL),
    data   = list(positions = as.data.frame(bpos)),
    class = "ipfEstimation"))

}

#' Estimates the location of the test observations
#'
#' @param ipfmodel  an ipfModel
#' @param test_fgp  a matrix or a data frame containing the fingerprints of the test set
#' @param test_pos  a matrix or a data frame containing the position of the test set fingerprints
#'
#' @return         An S3 object of class ipfEstimation, with the following properties:
#'                 location ->   a matrix with the predicted locations
#'                 errors   ->   a numeric vector with the errors
#'                 neighbors ->  a matrix with k columns and nrow(test) rows, with the
#'                               k most similar training observation for each test observation
#'                 weights ->    a matrix with k columns and nrow(test) rows, with the weights
#'
#' @examples
#'
#'     model <- ipfKnn(ipftrain[1:200, 1:168], ipftrain[1:200, 169:170])
#'     estimation <- ipfEstimate(model, ipftest[1:20, 1:168], ipftest[1:20, 169:170])
#'
#'\dontrun{
#'     model <- ipfProbabilistic(ipftrain[, 1:168], ipftrain[, 169:170], k = 9, delta = 10)
#'     estimation <- ipfEstimate(model, ipftest[, 1:170], ipftest[, 169:170])
#' }
#' @importFrom dplyr group_by arrange "%>%"
#'
#' @importFrom methods new
#'
#' @export
ipfEstimate <- function(ipfmodel, test_fgp, test_pos = NULL) {
  if (!inherits(ipfmodel, "ipfModel")) {
    stop("ipfEstimate: Wrong parameter type. ipfModel expected.")
  }

  extra_params <- ipfmodel$params$extra_params
  if (is.null(extra_params)) extra_params <- list()

  if (ipfmodel$params$name == 'prob') {
    est <- do.call(prob, c(list(
      tr_mns = ipfmodel$data$means,
      tr_sds = ipfmodel$data$sds,
      ts_fgp = test_fgp,
      k = ipfmodel$params$k,
      FUN = ipfmodel$params$FUN,
      delta = ipfmodel$params$delta
    ), extra_params))

  } else if (ipfmodel$params$name == 'knn') {
    est <- do.call(knn, c(list(
      tr_fgp = ipfmodel$data$fingerprints,
      ts_fgp = test_fgp,
      k = ipfmodel$params$k,
      method = ipfmodel$params$method,
      weights = ipfmodel$params$weights,
      norm = ipfmodel$params$norm,
      sd = ipfmodel$params$sd,
      epsilon = ipfmodel$params$epsilon,
      alpha = ipfmodel$params$alpha,
      threshold = ipfmodel$params$threshold,
      FUN = ipfmodel$params$FUN
    ), extra_params))
  } else if (ipfmodel$params$name == 'prox') {
    return (prox(
      fingerprints = test_fgp,
      fgps_pos = test_pos,
      bpos = ipfmodel$data$positions,
      rssirange = ipfmodel$params$rssirange,
      norssi = ipfmodel$params$norssi,
      alpha = ipfmodel$params$alpha,
      wapPow1 = ipfmodel$params$wapPow1
    ))
  } else {
    stop("Invalid model.")
  }

  tr_pos <- ipfmodel$data$positions
  if (!is.null(test_pos)) {
    if (!is.data.frame(test_pos)) {
      test_pos <- data.frame(test_pos)
    }
    test_pos <- data.frame(test_pos[1:nrow(test_fgp),])
  }

  col_errors <- rep(0, nrow(test_fgp))
  mloc <- matrix(0, nrow(test_fgp), ncol(tr_pos))
  mloc <- data.frame(mloc)

  n_factors <- sum(sapply(tr_pos, function(x) is.factor(x)))
  confusion <- vector('list', n_factors)
  nms <- vector('character', n_factors)
  nc <- 1

  for (i in 1:ncol(mloc)) {
    if (is.factor(tr_pos[, i])) {
      fm <- matrix(tr_pos[est$neighbors, i], ncol = ipfmodel$params$k)
      mloc[, i] <- factor(apply(fm, 1, function(x) names(which.max(table(x)))))
      if (!is.null(test_pos)) {
        name <- names(tr_pos)[i]
        if (is.null(name)) {
          name = paste0('confusion', nc)
        }
        confusion[[nc]] <- table(test_pos[, i], mloc[, i])
        nms[nc] <- name
        nc <- nc + 1
      }
    } else {
      mloc[, i] <- rowSums(tr_pos[est$neighbors, i] * est$weights)
      if (!is.null(test_pos)) {
        col_errors <- col_errors + (mloc[, i] - test_pos[, i]) ^ 2
      }
    }
  }

  colnames(mloc) <- colnames(tr_pos)
  names(confusion) <- nms

  errors <- c()
  if (!is.null(test_pos)) {
    errors <- sqrt(col_errors)
  }

  return(structure(list(
    location = mloc,
    errors = errors,
    confusion = confusion,
    neighbors = est$neighbors,
    weights = est$weights),
    class = "ipfEstimation"))
}


#' @importFrom stats optim
minim <- function(par, fingerprint, wapPos, wP1, alpha) {

  # Escoger un punto inicial. El correspondiente al WAP con mayor señal.
  par = wapPos[which.max(fingerprint),]

  d1 <- 10 ^ ((wP1 - fingerprint) / (10 * alpha))
  sdd <- 4
  sd <- (sdd * log(10) / (10 * alpha))
  w <- 1 / (d1^2 * exp(sd ^ 2) * (exp(sd ^ 2) - 1))
  optim(par = par, method = 'BFGS', fn = df, d1 = d1, wapPos = wapPos, w = w)$par
}

# Function to minimize to find position from fingerprint and WAPs location
df <- function(pos, d1, wapPos, w) {
  d2 <- apply(wapPos, 1, function(x) sqrt((x[1] - pos[1])^2 + (x[2] - pos[2])^2))
  sum(unlist(sapply(1:ncol(d1), function(i) w[i] * (d1[,i] - d2[i])^2)), na.rm = TRUE)
}


# AVERAGE NUMBER OF WAPS PER OBSERVATION
rmce_getANW <- function(traindata, noRSSI = NA) {
  if (is.na(noRSSI)) {
    traindata[!is.na(traindata)] = 1
  } else {
    traindata[traindata != noRSSI] = 1
  }
  return(mean(rowSums(traindata, na.rm = TRUE)))
}

# STANDARD DEVIATION OF THE NUMBER OF WAPS PER OBSERVATION
rmce_getSDNW <- function(traindata, noRSSI = NA) {
  if (is.na(noRSSI)) {
    traindata[!is.na(traindata)] = 1
  } else {
    traindata[traindata != noRSSI] = 1
  }
  return(sd(rowSums(traindata, na.rm = TRUE)))
}

# AVERAGE RSSI (WHEN > 0)
rmce_getAIW <- function(traindata, noRSSI = NA) {
  if (is.na(noRSSI)) {
    return(mean(c(t(traindata[!is.na(traindata)]))))
  } else {
    return(mean(c(t(traindata[traindata != noRSSI]))))
  }
}

# MEDIAN RSSI (WHEN > 0)
rmce_getMIW <- function(traindata, noRSSI = NA) {
  if (is.na(noRSSI)) {
    return(median(c(t(traindata[!is.na(traindata)]))))
  } else {
    return(median(c(t(traindata[traindata != noRSSI]))))
  }
}

# STANDARD DEVIATION OF RSSI (WHEN > 0)
rmce_getSDIW <- function(traindata, noRSSI = NA) {
  if (is.na(noRSSI)) {
    return(sd(c(t(traindata[!is.na(traindata)]))))
  } else {
    return(sd(c(t(traindata[traindata != noRSSI]))))
  }
}

# NUMBER OF OBSERVATIONS
rmce_getNP <- function(locdata) {
  return(nrow(locdata))
}

# NUMBER OF GROUPS
rmce_getNG <- function(locdata) {
  return(length(unique(ipfGroup(locdata))))
}

# AVERAGE DENSITY OF OBSERVATIONS
rmce_getADP <- function(locdata) {
  return(mean(rowMeans(ipfDistance(locdata, locdata), na.rm = TRUE)))
}

# AVERAGE DENSITY OF GROUPS
#' @importFrom dplyr group_by arrange "%>%"
#'
rmce_getADG <- function(locdata) {
  meanNA <- function(x) mean(x, na.rm = TRUE)
  GROUP <- NULL
  gr <- ipfGroup(locdata)
  locdata$GROUP <- gr
  locdata <- locdata %>%
    group_by(GROUP) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), meanNA), .groups = "drop") %>%
    arrange(GROUP)
  locdata$GROUP <- NULL
  locdata <- as.data.frame(locdata)
  return(rmce_getADP(locdata))
}

# AVERAGE DISTANCE BETWEEN TWO RANDOM POINTS INSIDE A RECTANGLE
rmce_getAAD <- function(locdata) {
  return(rectAD(rmce_getW(locdata), rmce_getH(locdata)))
}

# WIDTH
rmce_getW <- function(locdata) {
  minW <- min(locdata[, 1], na.rm = TRUE)
  maxW <- max(locdata[, 1], na.rm = TRUE)
  return(abs(maxW - minW))
}

# HEIGHT
rmce_getH <- function(locdata) {
  minH <- min(locdata[, 2], na.rm = TRUE)
  maxH <- max(locdata[, 2], na.rm = TRUE)
  return(abs(maxH - minH))
}

# TRAINING RADIO MAP DATA PER CELL
rm_getCD <- function(traindata, locdata, noRSSI = NA, gridSize = 5) {
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

      nps <- c(nps, nrow(p))                   # Number of observations in the grid
      ngs <- c(ngs, rmce_getNG(p))             # Number of groups in the grid
      anw <- c(anw, rmce_getANW(t, noRSSI))    # Average number of WAPs heard per observation in the grid
      aiw <- c(aiw, rmce_getAIW(t, noRSSI))    # Average RSSI in the grid
      miw <- c(miw, rmce_getMIW(t, noRSSI))    # Median RSSI in the grid
      sdw <- c(sdw, rmce_getSDIW(t, noRSSI))   # Stardard deviation of RSSI in the grid
    }
  }

  return(list(
    anppc = mean(nps), angpc = mean(ngs), sdnppc = sd(nps), sdngpc = sd(ngs),
    anwpc = mean(anw), sdnwpc = sd(anw), aiwpc = mean(aiw), sdiwpc = sd(aiw),
    area = ng * gridSize * gridSize
  ))
}


rmce_getTCV <- function(fingerprints, positions, rangeRSSI, noRSSI) {
  set.seed(42)
  fingerprints <- ipfTransform(fingerprints, inRange = rangeRSSI, outRange = c(0, 100),
                               inNoRSSI = noRSSI, outNoRSSI = 0)
  N       <- 10
  cv      <- 0.7
  error   <- 0
  tot_n   <- nrow(fingerprints)
  tr_size <- floor(tot_n * cv)
  for (i in 1:N) {
    tr_ind   <- sample(1:tot_n, tr_size)

    tr_set_f <- fingerprints[tr_ind, ]
    ts_set_f <- fingerprints[-tr_ind, ]
    tr_set_p <- positions[tr_ind, ]
    ts_set_p <- positions[-tr_ind, ]

    knnModel <- ipfKnn(tr_set_f, tr_set_p, k = 1)
    knnEst   <- ipfEstimate(knnModel, ts_set_f, ts_set_p)
    error    <- error + mean(knnEst$errors)
  }

  result <- error / N
  return(result)
}


ipfrmdata <- function(fingerprints, positions, rangeRSSI, noRSSI, gridSize) {
  ANW  <- rmce_getANW(fingerprints, noRSSI)                       # AVERAGE NUMBER OF WAPS PER OBS
  AIW  <- rmce_getAIW(fingerprints, noRSSI)                       # AVERAGE RSSI
  MIW  <- rmce_getMIW(fingerprints, noRSSI)                       # MEDIAN RSSI
  SDIW <- rmce_getSDIW(fingerprints, noRSSI)                      # STANDARD DEVIATION OF RSSI
  SDNW <- rmce_getSDNW(fingerprints, noRSSI)                      # STANDARD DEVIATION OF RSSI
  TCV  <- rmce_getTCV(fingerprints, positions, rangeRSSI, noRSSI) # STANDARD DEVIATION OF RSSI

  NP   <- rmce_getNP(positions)    # TOTAL NUMBER OF OBSERVATIONS
  NG   <- rmce_getNG(positions)    # TOTAL NUMBER OF GROUPS
  ADP  <- rmce_getADP(positions)   # AVERAGE DISTANCE BETWEEN OBSERVATIONS
  ADG  <- rmce_getADG(positions)   # AVERAGE DISTANCE BETWEEN GROUPS
  AAD  <- rmce_getAAD(positions)   # AVERAGE DISTANCE BETWEEN TWO RANDOM POINTS (RECTANGLE WxH)
  W    <- rmce_getW(positions)     # TOTAL AREA WIDTH
  H    <- rmce_getH(positions)     # TOTAL AREA HEIGHT


  rmcd <- rm_getCD(fingerprints, positions, noRSSI = noRSSI, gridSize = GRIDSIZE)

  return(structure(list(
    ANW = ANW,
    AIW = AIW,
    MIW = MIW,
    SDIW = SDIW,
    SDNW = SDNW,
    TCV = TCV,
    NP = NP,
    NG = NG,
    ADP = ADP,
    ADG = ADG,
    AAD = AAD,
    W = W,
    H = H,
    anppc = rmcd$anppc,
    angpc = rmcd$angpc,
    sdnppc = rmcd$sdnppc,
    sdngpc = rmcd$sdngpc,
    anwpc = rmcd$anwpc,
    sdnwpc = rmcd$sdnwpc,
    aiwpc = rmcd$aiwpc,
    sdiwpc = rmcd$sdiwpc,
    area = rmcd$area
  ), class = "rmdata"))

}



#' Estimates the inherent difficulty of the radio map
#'
#' @param fingerprints  a matrix or a data frame containing the RSSI data (fingerprints) of the
#'                      observations
#' @param positions     a matrix or a data frame containing the positions of the fingerprints
#' @param rangeRSSI     range of the RSSI data
#' @param noRSSI        value used to represent a not detected AP
#' @param gridSize      size of the grid to consider
#'
#' @return         a numeric value representing the RMID value (Radio Map Inherent Difficulty)
#'
#' @examples
#'\dontrun{
#'     rmid <- ipfRMID(ipftrain[1:200, 1:168], ipftrain[1:200, 169:170], noRSSI = NA)
#'}
#' @export
ipfRMID <- function(fingerprints, positions, rangeRSSI = c(-100, 0), noRSSI = NA,
                    gridSize = 5) {

  rmdata <- ipfrmdata(fingerprints, positions, rangeRSSI, noRSSI, gridSize)

  b <- rmdata$ADP * rmdata$SDNW * rmdata$SDIW * rmdata$TCV
  a <- 10000
  rmid <- b / (a + b)
  return(rmid)
}



rectAD <- function(L1, L2) {
  d <- sqrt(L1^2 + L2^2)
  t1 <- (L1^3) / (L2^2) + (L2^3) / (L1^2)
  t2 <- d * (3 - (L1^2) / (L2^2) - (L2^2) / (L1^2))
  t3 <- (5 / 2) * (((L2^2) / L1) * log((L1 + d) / L2) + ((L1^2) / L2) * log((L2 + d) / L1))
  return((t1 + t2 + t3) / 15)
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
#'     model <- ipfKnn(ipftrain[1:200, 1:168], ipftrain[1:200, 169:170])
#'     estimation <- ipfEstimate(model, ipftest[1:20, 1:168], ipftest[1:20, 169:170])
#'     ipfPlotPdf(estimation)
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density geom_vline scale_color_manual
#' @importFrom ggplot2 element_rect scale_alpha_manual scale_x_continuous labs theme element_text
#' @importFrom methods is
#'
#' @export
ipfPlotPdf <- function(estimation, xlab = 'error', ylab = 'density',
                       title = 'Probability density function') {
  if (inherits(estimation, "ipfEstimation")) {
    if (is.null(estimation$errors)) {
      stop("The estimation has no errors data.")
    }
  } else {
    stop("ipfPlotPdf: Wrong parameter type. ipfEstimation expected.")
  }

  # Avoid no visible binding for global variable
  ..density.. <- NULL
  median <- NULL
  errors <- data.frame(ERRORS = estimation$errors)

  tickroundup <- c(1, 2, 5, 10)
  xmax <- max(errors)
  xmark <- xmax / 20
  xtick <- 10^floor(log10(xmark)) *
    tickroundup[[which(xmark <= 10^floor(log10(xmark)) * tickroundup)[[1]]]]

  ggplot(data = errors, aes(errors$ERRORS)) +
    geom_histogram(aes(y = ..density..), binwidth = max(errors$ERRORS) / 50, alpha = .5) +
    geom_density(alpha = .2, fill = "#666666", colour = 'grey30') +
    geom_vline(aes(xintercept = median(errors$ERRORS), color = "Median", linetype = "Median"),
               alpha = 0.5, size = 1, show.legend = TRUE) +
    geom_vline(aes(xintercept = mean(errors$ERRORS), color = "Mean", linetype = "Mean"),
               alpha = 0.5, size = 1, show.legend = TRUE) +
    scale_color_manual(name = '', values = c(Median = "royalblue1", Mean = "firebrick1")) +
    scale_alpha_manual(name = '', values = c(0.5, 0.5)) +
    scale_x_continuous(breaks = seq(0, xmax, xtick)) +
    scale_linetype_manual(name = '', values = c(Median = "longdash", Mean = "dotted")) +
    labs(x = "error",y = "density")  +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.9, 0.85),
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
#'     model <- ipfKnn(ipftrain[1:200, 1:168], ipftrain[1:200, 169:170])
#'     estimation <- ipfEstimate(model, ipftest[1:20, 1:168], ipftest[1:20, 169:170])
#'     ipfPlotEcdf(estimation)
#'
#' @importFrom ggplot2 ggplot aes stat_ecdf geom_vline scale_color_manual element_text
#' @importFrom ggplot2 element_rect scale_alpha_manual scale_x_continuous labs theme scale_linetype_manual
#' @importFrom methods is
#' @importFrom stats median pnorm sd
#'
#' @export
ipfPlotEcdf <- function(estimation, xlab = 'error',
                        ylab = 'cumulative density of error',
                        title = 'Empirical cumulative density function') {
  if (inherits(estimation, "ipfEstimation")) {
    if (is.null(estimation$errors)) {
      stop("The estimation has no errors data.")
    }
  } else {
    stop("ipfPlotEcdf: Wrong parameter type. ipfEstimation expected.")
  }

  errors <- data.frame(ERRORS = estimation$errors)

  tickroundup <- c(1, 2, 5, 10)
  xmax <- max(errors)
  xmark <- xmax / 20
  xtick <- 10^floor(log10(xmark)) *
    tickroundup[[which(xmark <= 10^floor(log10(xmark)) * tickroundup)[[1]]]]

  ggplot(errors, aes(errors$ERRORS)) +
    stat_ecdf(geom = "line", alpha = .75) +
    geom_vline(aes(xintercept = median(errors$ERRORS), color = "Median", linetype = "Median"),
               alpha = 0.5, size = 1, show.legend = TRUE) +
    geom_vline(aes(xintercept = mean(errors$ERRORS), color = "Mean", linetype = "Mean"),
               alpha = 0.5, size = 1, show.legend = TRUE) +
    scale_color_manual(name = '', values = c(Median = "royalblue1", Mean = "firebrick1")) +
    scale_alpha_manual(name = '', values = c(0.5, 0.5)) +
    scale_x_continuous(breaks = seq(0, xmax, xtick)) +
    scale_linetype_manual(name = '', values = c(Median = "longdash", Mean = "dotted")) +
    labs(y = "cumulative density of error", x = "error") +
    labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = c(0.9, 0.85),
          legend.background = element_rect(color = "transparent",
                                           fill = "transparent",
                                           size = 0,
                                           linetype = "blank"))

}

#' Plots the spatial location of the observations
#'
#' @param positions   a data frame or matrix with the positions
#' @param plabel      if TRUE, adds labels to groups / observations
#' @param reverseAxis swaps axis
#' @param xlab        x-axis label
#' @param ylab        y-axis label
#' @param title       plot title
#' @param pgrid       plot grid (boolean)
#'
#' @examples
#'
#'     ipfPlotLocation(ipftrain[, 169:170])
#'
#'     ipfPlotLocation(ipftrain[, 169:170], plabel = TRUE, reverseAxis = TRUE,
#'                     title = 'Position of training set observations')
#'
#' @importFrom dplyr group_by arrange "%>%"
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_rect
#'
#' @export
ipfPlotLocation <- function(positions, plabel = FALSE, reverseAxis = FALSE, xlab = NULL, ylab = NULL,
                       title = '', pgrid = FALSE) {

  if (reverseAxis) {
    positions <- positions[,c(2, 1)]
  }

  p <- ggplot(positions)
  # To avoid no visible binding for global variable
  GROUP <- NULL
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL

  if (pgrid) {
    minX <- min(positions[, 1])
    maxX <- max(positions[, 1])
    minY <- min(positions[, 2])
    maxY <- max(positions[, 2])

    columns <- ceiling(abs(maxX - minX) / GRIDSIZE)
    rows <- ceiling(abs(maxY - minY) / GRIDSIZE)
    rect.data <- NULL

    for (c in 1:columns) {
      for (r in 1:rows) {
        lxl <- minX + (c - 1) * GRIDSIZE
        lxr <- lxl + GRIDSIZE
        lyb <- minY + (r - 1) * GRIDSIZE
        lyt <- lyb + GRIDSIZE

        s <- positions[, 1] >= lxl & positions[, 1] < lxr &
             positions[, 2] >= lyb & positions[, 2] < lyt

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

  p <- p + geom_point(aes(x = positions[, 1], y = positions[, 2]), alpha = .25)

  if (plabel) {
    groups <- ipfGroup(positions)
    label.data <- positions
    label.data$GROUP = groups
    label.data <- data.frame(label.data %>%
      group_by(GROUP) %>%
      dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop") %>%
      arrange(GROUP))
    p <- p + geom_text(
      data = label.data,
      aes(x = label.data[, 2], y = label.data[, 3], label = label.data$GROUP),
      hjust = 0, vjust = 0, color = 'grey20', size = 3, nudge_x = 0.05, nudge_y = 0.05)
  }

  if (is.null(xlab)) {
    xlab <- colnames(positions)[1]
  }

  if (is.null(ylab)) {
    ylab <- colnames(positions)[2]
  }

  p <- p + labs(y = ylab, x = xlab) + labs(title = title)

  p

}

#' Plots the estimated locations
#'
#' @param model           an ipfModel
#' @param estimation      an ipfEstimation
#' @param testpos         position of the test observations
#' @param observations    a numeric vector with the indices of estimations to plot
#' @param reverseAxis     swaps axis
#' @param showneighbors   plot the k selected neighbors
#' @param showLabels      shows labels
#' @param xlab            x-axis label
#' @param ylab            y-axis label
#' @param title           plot title
#'
#' @examples
#'
#'     model      <- ipfKnn(ipftrain[1:200, 1:168], ipftrain[1:200, 169:170])
#'     estimation <- ipfEstimate(model, ipftest[1:20, 1:168], ipftest[1:20, 169:170])
#'     ipfPlotEstimation(model, estimation, ipftest[1:20, 169:170],
#'                       observations = 1:4, showneighbors = TRUE,
#'                       reverseAxis = TRUE)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_curve arrow unit
#'
#' @export
ipfPlotEstimation <- function(model, estimation, testpos = NULL, observations = c(1),
                              reverseAxis = FALSE, showneighbors = FALSE, showLabels = FALSE,
                              xlab = NULL, ylab = NULL, title = '') {

  ePoints <- as.data.frame(matrix(as.matrix(estimation$location[observations, ]), length(observations), 2))
  if (reverseAxis) {
    if (!is.null(testpos)) {
      testpos <- testpos[,c(2, 1)]
    }
    ePoints <- ePoints[,c(2, 1)]
  }

  sPoints <- NULL

  p <- ggplot() + geom_point(data = ePoints, aes(x = ePoints[, 1], y = ePoints[, 2]),
                             size = 6, alpha = 1, shape = 0)

  if (showneighbors) {
    nPoints <- NULL
    for (i in 1:length(observations)) {
      np <- model$data$positions[estimation$neighbors[observations[i],], 1:2]
      if (reverseAxis) {
        np <- np[,c(2, 1)]
      }
      nPoints <- rbind(nPoints, np)
      for (j in 1:model$params$k) {
        sPoints <- rbind(sPoints, cbind(ePoints[i, ], nPoints[j + (i - 1) * model$params$k, ]))
      }
    }

    p <- p + geom_point(data = nPoints, aes(x = nPoints[, 1], y = nPoints[, 2]),
                        size = 5, col = 'black', shape = 4) +
             geom_segment(data = sPoints, aes(x = sPoints[,1], xend = sPoints[, 3],
                                              y = sPoints[, 2], yend = sPoints[, 4]),
               alpha = 0.5, col = 'blue', size = 1, linetype = 'dashed')
  }


  if (!is.null(testpos)) {
    cPoints <- NULL
    for (i in 1:length(observations)) {
      cPoints <- rbind(cPoints, cbind(ePoints[i, ], testpos[observations[i], ]))
    }
    cPoints[cPoints[,1] == cPoints[,3] & cPoints[,2] == cPoints[,4], ] <- NA
    p <- p + geom_point(data = cPoints,
                       aes(x = cPoints[, 3], y = cPoints[, 4]),
                       size = 4, alpha = 1, col = 'blue', shape = 1,  na.rm = TRUE)
    p <- p + geom_curve(data = cPoints,
                        aes(x = cPoints[,3], xend = cPoints[, 1],
                            y = cPoints[, 4], yend = cPoints[, 2]),
                        alpha = 0.7, col = 'red', size = 1,
                        arrow = arrow(length = unit(0.02, "npc")),
                        curvature = 0.2, na.rm = TRUE)
  }

  if (showLabels) {
    ePoints$LABELS <- observations
    p <- p + geom_text(
      data = ePoints,
      aes(x = ePoints[, 1], y = ePoints[, 2], label = ePoints$LABELS),
      hjust = -0.5, vjust = -0.5, color = 'black', size = 4)
  }

  if (is.null(xlab)) {
    if (!is.null(testpos)) {
      xlab <- colnames(testpos)[1]
    } else {
      xlab <- ''
    }
  }

  if (is.null(ylab)) {
    if (!is.null(testpos)) {
     ylab <- colnames(testpos)[2]
    } else {
      ylab <- ''
    }
  }
  p <- p + labs(y = ylab, x = xlab) + labs(title = title)

  return(p)
}



#' # Given a dataset and a distance, this function computes
#' # the clusters of the dataset, calculates the mean RSSI for
#' # all components of the cluster and returns a new dataset
#' # with a row (of means and centroid) per cluster
#' #' @importFrom leaderCluster leaderCluster
#' ipfClD <- function(locdata, distance) {
#'   cl <- leaderCluster(locdata, distance)
#'   clusters <- data.frame(row = row.names(locdata), clusterID = cl$cluster_id)
#'   return (list(clusters = clusters, centroids = cl$cluster_centroids))
#' }
