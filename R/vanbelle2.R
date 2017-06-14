#### Function to fit using the second version of ranking approach

# function vanbelle2Fit implements the first version of the ranking approach of svm for
# survival time analysis. This version is called ranksvmc in my master thesis
#
# @param X [matrix(1)]
#   Data set containing the training points
# @param Y [vector(1)]
#   vector of survival times
# @param delta [vector(1)]
#   statut vector: 1 for not censored
# @param kernel_type [character(1)]
#   indicates which kernel type will be used to build the kernel matrix
# @param kernel_pars [vector(1)]
#   vector containing the kernel parameter, when required. See the function kernel_matrix
#   for more details
# @param bin_cat [vector(1)]
#   set of index indicating which columns of the training set X muss treated as binar or
#   categorial varibales. Only required by additive kernel. See the function
#   kernel_matrix for more details
# @param makediff [function(1)]
#   names the function that muss be invoking to construct the matrix of differences on
#   all the comparable pairs
# @param opt_alg [character(1)]
#   tells which function muss be invoked to solve the quadratic optimisation problem
# @param  sgf_sv[integer(1)]
#   indicates how long the decimal part of the solutions muss be rounded
# @param sigf [integer(1)]
#   used by 'ipop' when required. See 'ipop' for more details
# @param maxiter [integer(1)]
#   used by 'ipop' when required. See 'ipop' for more details
# @param margin [integer(1)]
#   used by 'ipop' when required. See 'ipop' for more details
# @param bound [integer(1)]
#   used by 'ipop' when required. See 'ipop' for more details
# @return [VB2FitObj(1)]
#              alpha.fact [vector(1)] estimated factors
#              Xtrain [matrix(1)] number of support vectors
#              Dc [matrix(1)]
#              kernel_type [character] kernel used during fitting phase
#              kernel_pars [vector(1)] parameters used to construict the kernel matrix
#-----------------------------------------------------------------------------------------
#' fits the 'vanbelle2' version of the ranking approach of survival support vector ananlysis.
#'
#'
#' @title survivalsvm (ranking approach)
#' @param X [\code{matrix(1)}]\cr
#' Matrix of training data point.
#' @param Y [\code{vector(1)}]\cr
#' Vector of survival times.
#' @param delta [\code{vector(1)}]\cr
#' Vector of status: 1 = not censored.
#' @param meth_par [\code{numeric(1)}]\cr
#' Parameter of regularization.
#' @param kernel_type [\code{numeric(1)}]\cr
#' Kernel that will be used to fit the model. The handled type are: linear kern ('lin_kern'), additive kernel ('add_kernel'),
#' radial basis kernels ('rbf_kernel' and 'rbf4_kernel') and the polynomial kernel ('poly_kernel').
#' @param kernel_pars [\code{numeric(1)|vector(1)}]\cr
#' Parameters of kernel, when required.
#' @param bin_cat [\code{vector(1)}]\cr
#' Indexes of binary/categorical varibales
#' @param makediff [\code{character(1)}]\cr
#' String indicating which of \code{'makediff1'}, \code{'makediff2'} or \code{'makediff3'}
#' will be used.
#' @param opt_alg [\code{character}]\cr
#' Program that will be used to solve the quadratic optimization problem. Either \code{\link{quadprog}} or \code{\link{ipop}}.
#' @param sgf_sv [\code{integer(1)}]\cr
#' Number of decimal digits in the solution of the quadratic optimization problem.
#' @param sigf [\code{integer(1)}]\cr
#' Used by \code{\link{ipop}}. See \code{\link{ipop}} for details.
#' @param maxiter [\code{inter(1)}]\cr
#' Used by \code{\link{ipop}}. See \code{\link{ipop}} for details.
#' @param margin [\code{numeric(1)}]\cr
#' Used by \code{\link{ipop}}. See \code{\link{ipop}} for details.
#' @param bound [\code{numeric(1)}]\cr
#' Used by \code{\link{ipop}}. See \code{\link{ipop}} for details.
#' @param eig.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#' @param conv.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#' @param posd.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#'
#' @export
#' @return [\code{VB2FitObj(1)}]
#' Object of class \code{VB2FitObj} containing elements:
#' \tabular{ll}{
#'    \code{Alpha} \tab solution of the quadratic optimization problem, \cr
#'    \code{Xtrain} \tab matrix of training points,\cr
#'    \code{DifMat} \tab matrix used to maked differences between neighbor points, \cr
#'    \code{Kernel} \tab kernel matrix, an object of class \code{Kernel},\cr
#'    \code{OptMeth} \tab program used to solve the quadratic optimization problem.\cr
#'  }
#'
#' @author Cesaire J. K. Fouodo
#' @keywords internal
vanbelle2Fit <- function (X, Y, delta,
                         meth_par = 1, kernel_type = "lin_kernel",
                         kernel_pars = NA, bin_cat = integer(0),
                         makediff = makediff3, opt_alg = "quadprog", sgf_sv = 5,
                         sigf = 7, maxiter = 40, margin = 0.05, bound = 10,
                         eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08) {
  if (!(opt_alg %in% c("quadprog", "ipop"))) {
    stop("'opt_alg' must be either 'quadprog' or 'ipop'")
  }
  i.ord <- order(Y)
  Y <- Y[i.ord]
  delta <- delta[i.ord]
  X <- X[i.ord,]
  n <- length(Y)
  if (is.na(kernel_pars) & !(kernel_type == "add_kernel" || kernel_type == "lin_kernel")) {
    kernel_pars <- rep(1/ncol(X), ncol(X))
  }
  # Build the Kernel matrix
  Ker <- kernelMatrix(Xtrain = X, kernel_type = kernel_type, kernel_pars = kernel_pars, bin_cat = bin_cat)
  K <- getMat.Kernel(Ker)
  # Build the matrix D of comparables pairs
  md <- makediff(Y = Y, delta = delta)
  Dc <- getMat.Diffmatrix(md)
  D <- crossprod(t(Dc), tcrossprod(K, Dc))
  # Solves the dual problem
  opt <- if (opt_alg == "quadprog") {
    pracma::quadprog(C = D <- as.matrix(Matrix::nearPD(D, eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)$mat),
                     d = -drop(crossprod(t(Dc), Y)), A = -diag(n-1), b = rep(0, n-1), lb = 0, ub = meth_par)
  } else {
    kernlab::ipop(H =  D, c = -crossprod(t(Dc), Y), A = t(rep(1, n-1)), b = 0, l = rep(0, n-1), u = rep(meth_par, n-1),
         r = (n-1)*meth_par, sigf = sigf, maxiter = maxiter, margin = margin, bound = bound)
  }
  alphapar <- if(opt_alg == "ipop"){
    round(opt@primal, sgf_sv)
  } else {
    round(opt$xmin, sgf_sv)
  }
  vb2fo <- VB2FitObj(Alpha = alphapar, Xtrain = X, DifMat = md, Kernel = Ker, OptMeth = opt_alg)
  return(vb2fo)
}

#--- construct the VB2FitObj
#' Constructs object of class \code{VB2FitObj}.
#'
#' @title survivalsvm (ranking approach)
#' @param Alpha [\code{vector(1)}]\cr
#' Solution of the quadratic optimization problem of interest.
#' @param Xtrain [\code{matrix(1)}]\cr
#' Matrix of training data points.
#' @param DifMat [\code{\link{Diffmatrix}(1)}]\cr
#' Matrix used to maked differences between neighbor points.
#' @param Kernel [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#' @param OptMeth [\code{character(1)}]\cr
#' Program used to solve the optimization problem.
#'
#' @return \code{VB2FitObj}
#' Object of class \code{RegFitObj} containing elements:
#' \tabular{ll}{
#'    \code{Alpha} \tab solution of the quadratic optimization problem, \cr
#'    \code{Xtrain} \tab matrix of training points,\cr
#'    \code{DifMat} \tab matrix used to maked differences between neighbor points, \cr
#'    \code{Kernel} \tab kernel matrix, an object of class \code{Kernel},\cr
#'    \code{OptMeth} \tab program used to solve the quadratic optimization problem.\cr
#'  }
#' @keywords internal
VB2FitObj <- function (Alpha = NULL, Xtrain = NULL, DifMat = NULL, Kernel = NULL, OptMeth = NULL) {
  vb2o <- VB1FitObj(Alpha = Alpha, Xtrain = Xtrain, DifMat = DifMat, Kernel = Kernel, OptMeth = OptMeth)
  class(vb2o) <- append(class(vb2o), "VB2FitObj")
  return(vb2o)
}
