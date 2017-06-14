#' computes the concordance index.
#'
#' @title \code{cindex}
#' @param obj [\code{survivalsvmprediction}]\cr
#' Object of class \code{survivalsvmprediction}.
#' @param Y [\code{vector}(1)]\cr
#' A numeric vector of truth survival times obeserved.
#'
#' @export
#' @return [\code{Integer}]
#' Concordance index.
#' @keywords internal
conindex <- function(obj, Y){
  if (!inherits(obj, "survivalsvmprediction")) {
    stop("Error: 'obj' must an object af class survivalsvmprediction.")
  }
  if (is.null(obj$predicted)) {
    stop("Error: field 'predicted' not found in 'obj'.")
  }
  X <- obj$predicted
  if (length(X) != length(Y)) {
    stop("Error: lengths do not macht.")
  }
  ci = Hmisc::rcorr.cens(x = X, S = Y)
  return(ci["C Index"])
}

#' compute the Logrank
#'
#' @param t1 [\code{vector}(1)]\cr
#' A numeric vector.
#' @param d1 [\code{vector}(1)]\cr
#' Binary vector.
#' @param t2 [numeric(1)]\cr
#' A numeric vector.
#' @param d2 [\code{vector}(1)]\cr
#' A binary vector.
#' @export
#'
#' @return list of:
#' \tabular{ll}{
#'  \code{chi_sq} \tab chi-squared statistic at a significance level of 95 \% and one degree of freedom, \cr
#'  \code{chi_p} \tab chi-squared probality at a significance level of 95 \% and one degree of freedom. \cr
#' }
#' @keywords internal
logrank <- function(t1, d1, t2, d2){
  t1.ord <- order(t1)
  t1 <- t1[t1.ord]
  d1 <-  d1[t1.ord]
  t2.ord <- order(t2)
  t2 <- t2[t2.ord]
  d2 <-  d2[t2.ord]
  t <- c(t1, t2)
  d <- c(d1, d2)
  times <- unique(t[d == 1])
  n <- length(times)
  tel.ner <- sapply(times, function(i){
    o <- sum(t == i) # failures
    r <- sum(t >= i) # at risk
    o1 <- sum(t1 == i) # failures
    r1 <- sum(t1 >= i) # at risk
    o2 <- sum(t2 == i) # failures
    r2 <- sum(t2 >= i) # at risk
    teller <- o1 - r1*o/r
    noemer <- r2*r1*o*(r-o) / (r^2 *(r - 1))
    return(c(teller, noemer))
  })
  search.na <- colSums(tel.ner)
  na.index <- which(is.na(search.na))
  if(length(na.index) > 0){
    tel.ner <- tel.ner[, -na.index]
  }
  chi <- sum(tel.ner[1,])^2 / sum(tel.ner[2,])
  return(list(chi_sq = chi,
              chi_p = 1 - stats::pchisq(chi, 1)))
}

#' computes the logrank statistic.
#'
#' @title \code{getLogrank}
#' @param obj [\code{survivalsvmprediction}(1)]\cr
#' Object of class \code{survivalsvmprediction}.
#' @param t [\code{numeric}(1)]\cr
#' Numeric vector (of survival times).
#' @param delta [\code{vector}(1)]\cr
#' Binary vector (of status).
#' @export
#'
#' @return list of:
#' \tabular{ll}{
#'  \code{chi_sq} \tab chi-squared statistic and \cr
#'  \code{chi_p} \tab chi-squared probality. \cr
#' }
#' @keywords internal
getLogrank <- function(obj, t, delta){
  if (!inherits(obj, "survivalsvmprediction")) {
    stop("Error: 'obj' must be an object af class survivalsvmprediction.")
  }
  if (is.null(obj$predicted)) {
    stop("Error: field 'predicted' not found in 'obj'.")
  }
  u <- obj$predicted
  part2 <- which(u > mean(u, na.rm = TRUE))
  part1 <- setdiff(1:length(t), part2)
  return(logrank(t1 = t[part1], d1 = delta[part1], t2 = t[part2], d2 = delta[part2]))
}
