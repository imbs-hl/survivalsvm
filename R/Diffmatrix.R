# function makediff3 build the matrix form that will be use tu compute differences
#
#                    This is exactly the version explain in my master thesis
# @param Y [vector(1)]
#   Y is thevector of survival times
# @param delta[vector(1)]
#   delta is the vector containing the censoring information. 1: not censored.
# @return [matrix(1)]
#   matrix used to compute differences between comparable data points amongst observations (to make thing like: y_i - y_{\bar(j)(i)}).
#----------------------------------------------------------------------------------------
#'  The first observation is assumed to be not censored (delta = 1).
#'  The Difference is computed between data point \code{i}
#'  and its neighbour that has the largest survival time but smaller than \code{y_i}, the survival time of \code{i}.
#'
#' @title \code{Diffmatrix}
#' @param Y  [\code{vector(1)}]\cr
#' Ordered \code{vector} of survival times.
#' @param delta [\code{vector(1)}]\cr
#' Vector of status.
#'
#' @return [\code{\link{Diffmatrix}(1)}]
#' Object of class \code{Diffmatrix} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab type of makediff function used to compute differences between neighbours. \cr
#'    \code{Mat} \tab matrix used to compute differences between comparable data points. \cr
#'  }
#' @export
#'
#' @examples Y <- c(1,3,3.5,4,8); delta <- c(0,0,1,1,0); makediff3(Y, delta)
#' @seealso \code{\link{makediff1}} and \code{\link{makediff2}}
makediff3 <- function (Y, delta) {
  if (!prod(sort(Y) == Y)) {
    stop("sort Y first!")
    }
  delta[1] <- 1 # assuming that the first entry is an event
  n <- length(Y)
  D <- sapply(n:2, function(i){
    v <- numeric(n)
    v[i] <- 1
    small <- which(delta == 1 & Y <= Y[i])
    smallest <- max(setdiff(small, i:n))
    v[smallest] <- -1
    return(v)
  }, simplify = FALSE)
  D <- do.call("rbind", D)
  return(Diffmatrix(Type = "diff3", Mat = D[(n-1):1,]))
}

# function makediff2 build the matrix form that will be use tu compute differences
#
#                    This is another approach to build the matrix of differences, in
#                    which the comparaison only takes place between the not censored observations
# @param Y [vector(1)]
#   Is thevector of survival times.
# @param delta[vector(1)]
#   Is the vector containing the censoring information. 1: not censored.
# @return [matrix(1)]
#   Matrix of differences between observations (to make thing like: y_i - y_{\bar(j)(i)}).
#----------------------------------------------------------------------------------------
#'  computes the matrix difference only between not censored data points (delta = 1).
#'  The data points are asssumed to be sorted by survival time and the difference
#'  is computed only if both the comparable data points are not censored.
#'
#' @title \code{Diffmatrix}
#' @param Y  [\code{vector(1)}]\cr
#' Ordered \code{vector} of survival times.
#' @param delta [\code{vector(1)}]\cr
#' Vector of status
#'
#' @return [\code{\link{Diffmatrix}(1)}]
#' Object of class \code{Diffmatrix} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab type of makediff function used to compute differences between neighbours. \cr
#'    \code{Mat} \tab matrix used to compute differences between comparable data points. \cr
#'  }
#' @export
#'
#' @examples Y <- c(1,3,3.5,4,8); delta <- c(0,0,1,1,0); makediff2(Y, delta)
#' @seealso \code{\link{makediff1}} and \code{\link{makediff3}}
makediff2 <- function (Y, delta) {
  # supposes all entries to be no-censored when no delta's parameter
  n <- length(Y)
  if (is.unsorted(Y)) {
    stop("sort before invoking makeDiff")
  }
  events <- which(delta == 1)
  if (length(events) == 0) return(NULL)
  comp.i <- expand.grid(events, events)
  comp.i <- unique(t(apply(comp.i, 1, sort)))
  nr <- nrow(comp.i)
  no.const <- sapply(1:nr, function(i) {
    if (!length(unique(comp.i[i,])) == 1) {
      return(i)
    }
  })
  no.const <- no.const[lengths(no.const) > 0]
  comp.i <- comp.i[unlist(no.const), , drop = FALSE]
  nr <- nrow(comp.i)
  comp.m <- apply(comp.i, 1, function(i) {
    v <- numeric(n)
    v[i] <- c(-1, 1)
    return(v)
  })
  return(Diffmatrix(Type = "diff2", Mat = t(comp.m)))
}

# function makediff1 build the matrix form that will be use tu compute differences
#
#                    This is another approach to build the matrix of differences, in
#                    which the comparaison only takes place between two consecutive
#                    (order is the survival time) observations when the first one is not
#                    censored. Two consecutive censored observations are not comparable.
# @param Y [vector(1)]
#   Is thevector of survival times
# @param delta[vector(1)]
#   Ss the vector containing the censoring information. 1: not censored
# @return [matrix(1)]
#   Matrix of differences between observations (to make thing like: y_i - y_{\bar(j)(i)})
# ------------------------------------------------------------------------------------------------
#'  The data points are asssumed to be sorted by survival time. The comparison only takes place between two consecutivee observations when the first one is not
#'  censored (delta = 1).
#'
#' @title \code{Diffmatrix}
#' @param Y  [\code{vector(1)}]\cr
#' Ordered \code{vector} of survival times.
#' @param delta [\code{vector(1)}]\cr
#' Vector of status
#'
#' @return [\code{\link{Diffmatrix}(1)}]
#' Object of class \code{Diffmatrix} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab type of makediff function used to compute differences between neighbours. \cr
#'    \code{Mat} \tab matrix used to compute differences between comparable data points. \cr
#'  }
#' @export
#'
#' @examples Y <- c(1,3,3.5,4,8); delta <- c(0,0,1,1,0); makediff1(Y, delta)
#' @seealso \code{\link{makediff2}} and \code{\link{makediff3}}
makediff1 <- function (Y, delta) {
  n <- length(Y)
  D1 <- diag(n)
  D2 <- diag(n)
  D <- D1[2:n,] - D2[1:(n-1),]
  I <- (D %*% Y > 0) & (delta[1:(n-1)] == 1)
  R <- diag(length(I))
  diag(R) <- I
  Dc <- R %*% D
  return(Diffmatrix(Type = "diff1", Mat = Dc))
}

# >>>> contuction of Diffmatrix class
#--- the constructor
# ----------------------------------------------------------------------------
#' constructs objects of class \code{Diffmatrix}.
#'
#'
#' @title \code{Diffmatrix}.
#' @param Type  [\code{character(1)}]\cr
#' Indicates which difference is performed. This must be one of \code{\link{makediff1}}, \code{\link{makediff2}} other \code{\link{makediff3}}.
#' @param Mat [\code{matrix(1)}]\cr
#' Matrix used to perfom differences.
#'
#' @return [\code{Diffmatrix(1)}]
#' Mutated object of class \code{Diffmatrix} containing elements:\cr
#'  \tabular{ll}{
#'    \code{Type} \tab type of differences bildet between neighbors and \cr
#'    \code{Mat} \tab matrix used to perform differences between comparable data points. \cr
#'  }
#' @keywords internal
Diffmatrix <- function(Type = NULL, Mat = NULL) {
  d <- list(Type = Type, Mat = Mat)
  class(d) <- "Diffmatrix"
  return(d)
}

#--- Specific mutator method for the Diffmatrix type
#'  Mutator for the type of objects of class \code{Diffmatrix}
#'
#'
#' @title \code{Diffmatrix}
#' @param dm [\code{\link{Diffmatrix}(1)}]\cr
#' Object of class \code{Diffmatrix}.
#' @param dmtype [\code{character(1)}]\cr
#'  New type
#'
#' @return [\code{\link{Diffmatrix}(1)}]
#' Mutated object of class \code{Diffmatrix} containing elements.
#'
#' @keywords internal
#'
setType.Diffmatrix <- function(dm, dmtype) {
  dm$Type <- dmtype
  return(dm)
}

#--- Specific mutator method for the \code{Diffmatrix} matrix
#'  Mutator for the field \code{Mat} of objects of class \code{Diffmatrix}.
#'
#'
#' @title \code{Diffmatrix}.
#' @param dm [\code{\link{Diffmatrix}(1)}]\cr
#' Object of class \code{Diffmatrix}.
#' @param dmmat  [\code{matrix(1)}]\cr
#' Matrix used to perform differences between comparable data points.
#' @return [\code{\link{Diffmatrix}(1)}]
#' Mutated object of class \code{Diffmatrix} containing elements:
#'
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setMat.Diffmatrix <- function(dm, dmmat) {
  dm$Mat <- dmmat
  return(dm)
}

#>>>> Mutators for the kernel object

#' Mutator for objects of class \code{Diffmatrix}.
#'
#'
#' @title \code{Diffmatrix}.
#' @param dm [\code{\link{Diffmatrix}(1)}]\cr
#' Object of class \code{Diffmatrix}.
#'
#' @return [\code{character(1)}]
#' Type of Diffmatrix
#' @keywords internal
getType.Diffmatrix <- function(dm) {
  return(dm$Type)
}

#' To get the matrix used to perform differences between comparable data points for object of class Diffmatrix.
#'
#'
#' @title \code{Diffmatrix}.
#' @param dm [\code{\link{Diffmatrix}(1)}]\cr
#' Object of class \code{Diffmatrix}.
#'
#' @return The matrix used to perform differences between comparable data points.
#' @keywords internal
getMat.Diffmatrix <- function(dm) {
  return(dm$Mat)
}
