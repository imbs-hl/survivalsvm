#### Function to fit using the hybrid approach

# function hybridFit implements hybrid approach of svm for
# survival time analysis. This version is called ranksvmc in my master thesis
#
# @param X [matrix(1)]
#   Data set containing the training points
# @param Y [vector(1)]
#   Vector of survival times
# @param delta [vector(1)]
#   Statut vector: 1 for not censored
# @param kernel_type [character(1)]
#   Indicates which kernel type will be used to build the kernel matrix
# @param kernel_pars [vector(1)]
#   Vector containing the kernel parameter, when required. See the function kernel_matrix
#   for more details
# @param bin_cat [vector(1)]
#   Set of index indicating which columns of the training set X muss treated as binar or
#   categorial varibales. Only required by additive kernel. See the function
#   kernel_matrix for more details
# @param makediff [function(1)]
#   Names the function that muss be invoking to construct the matrix of differences on
#   all the comparable pairs
# @param opt_alg [character(1)]
#   Tells which function muss be invoked to solve the quadratic optimisation problem
# @param  sgf_sv[integer(1)]
#   Indicates how long the decimal part of the solutions muss be rounded
# @param sigf [integer(1)]
#   Used by 'ipop' when required. See 'ipop' for more details
# @param maxiter [integer(1)]
#   Used by 'ipop' when required. See 'ipop' for more details
# @param margin [integer(1)]
#   Used by 'ipop' when required. See 'ipop' for more details
# @param bound [integer(1)]
#   Used by 'ipop' when required. See 'ipop' for more details
# @return [VB2FitObj(1)]
#              alpha.fact [vector(1)] estimated factors
#              Xtrain [matrix(1)] number of support vectors
#              Dc [matrix(1)]
#              kernel_type [character] kernel used during fitting phase
#              kernel_pars [vector(1)] parameters used to construict the kernel matrix

#### Hybrid approach: the function fits the model2 #######
#-----------------------------------------------------------------------------------------
#' fits survivalsvm model based on hybrid approach method for survival support vector ananlysis.
#'
#'
#' @title survivalsvm (hybrid approach)
#' @param X [\code{matrix(1)}]\cr
#' Matrix of training data points.
#' @param Y [\code{vector(1)}]\cr
#' Vector of survival times.
#' @param delta [\code{vector(1)}]\cr
#' Vector of status: 1 = not censored.
#' @param meth_par [\code{numeric(1)}]\cr
#' Parameters of regularization.
#' @param kernel_type [\code{character(1)}]\cr
#' Kernel that is used to fit the model. The handled type are: linear kern ('lin_kern'), additive kernel ('add_kernel'),
#' radial basis kernels ('rbf_kernel' and 'rbf4_kernel') and the polynomial kernel ('poly_kernel').
#' @param kernel_pars [\code{numeric(1)|vector(1)}]\cr
#' Parameters of kernel, when required.
#' @param bin_cat [\code{vector(1)}]\cr
#' Indexes of binary/categorical varibales
#' @param makediff [\code{character(1)}]\cr
#' String indicating which of \code{'makediff1'}, \code{'makediff2'} or \code{'makediff3'}
#' will be used.
#' @param opt_alg [\code{character(1)}]\cr
#' Program that will be used to solve the quadratic optimization problem. Either \code{\link[pracma]{quadprog}} or \code{\link[kernlab]{ipop}}.
#' @param sgf_sv [\code{integer(1)}]\cr
#' Number of decimal digits in the solution of the quadratic optimization problem.
#' @param sigf [\code{integer(1)}]\cr
#' Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.
#' @param maxiter [\code{integer(1)}]\cr
#' Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.
#' @param margin [\code{numeric(1)}]\cr
#' Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.
#' @param bound [\code{numeric(1)}]\cr
#' Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.
#' @param eig.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.
#' @param conv.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.
#' @param posd.tol [\code{numeric(1)}]\cr
#' Used by \code{nearPD} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.
#'
#' @export
#' @return [\code{Hybrid(1)}]
#' Object of class \code{Hybrid} containing elements:
#' \tabular{ll}{
#'    \code{Alpha} \tab Solution of the quadratic optimization problem, \cr
#'    \code{Xtrain} \tab Matrix of training points,\cr
#'    \code{DifMat} \tab Matrix used to maked differences between neighbor points, \cr
#'    \code{Kernel} \tab Kernel matrix, an object of class \code{Kernel},\cr
#'    \code{OptMeth} \tab Program used to solve the quadratic optimization problem.\cr
#'  }
#'
#' @importFrom pracma quadprog
#' @importFrom kernlab ipop
#' @keywords internal
#' @author Cesaire J. K. Fouodo
hybridFit <- function (X, Y, delta,
                      meth_par = c(1, 1), kernel_type = "lin_kernel",
                      kernel_pars = NA, bin_cat = integer(0),
                      makediff = makediff3, opt_alg = "quadprog", sgf_sv = 5,
                      sigf = 7, maxiter = 40, margin = 0.05, bound = 10,
                      eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08) {
  # Implementation of the hybrid approach
  # Use censored failure time data
  # event(not censored): delta=1;  censored: delta=0
  if (!(opt_alg %in% c("quadprog", "ipop"))) {
    stop("opt_alg must be either 'quadprog' or 'ipop'")
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
  # Solve dual problem
  # Let first construct the Hessian matrix
  dm <- makediff(Y = Y, delta = delta)
  Dc <- getMat.Diffmatrix(dm)
  KDc.t <- tcrossprod(K, Dc)
  bloc1 <- crossprod(t(Dc), KDc.t)
  delta.KDc.t <- tcrossprod(delta * K, Dc)
  l1 <- cbind(crossprod(t(Dc), KDc.t), t(KDc.t), - t(delta.KDc.t))
  l2 <- cbind(KDc.t, K, -K)
  l3 <- cbind(-delta.KDc.t, -K, K)
  H <- rbind(l1, l2, l3)
  opt <- if (opt_alg == "quadprog") {
    pracma::quadprog(C =  as.matrix(Matrix::nearPD(H, eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)$mat),
                     d = c(-crossprod(t(Dc), Y), -Y, delta*Y),
             A = -diag(3*n - 1), b = rep(0, 3*n-1), lb = rep(0, 3*n-1), ub = c(rep(meth_par[1], n-1), rep(meth_par[2], 2*n)),
             Aeq = c(rep(0, n-1), rep(-1, n), delta), beq = 0)
  } else {
    kernlab::ipop(H = H, c = matrix(c(-crossprod(t(Dc), Y), -Y, delta*Y)), A = t(c(rep(0, n-1), rep(-1, n), delta)), b = 0,
         r = 0, l = matrix(rep(0, 3*n-1)), u = matrix(c(rep(meth_par[1], n-1), rep(meth_par[2], 2*n))), sigf = sigf,
         maxiter = maxiter, margin = margin, bound = bound)
  }
  if(FALSE){
    opt1 <- pracma::quadprog(C =  as.matrix(nearPD(H)$mat), d = c(-crossprod(t(Dc), Y), -Y, delta*Y),
                     A = -diag(3*n - 1), b = rep(0, 3*n-1), lb = rep(0, 3*n-1), ub = c(rep(meth_par[1], n-1), rep(meth_par[2], 2*n)),
                     Aeq = c(rep(0, n-1), rep(-1, n), delta), beq = 0)
    opt2 <- kernlab::ipop(H = H, c = matrix(c(-crossprod(t(Dc), Y), -Y, delta*Y)), A = t(c(rep(0, n-1), rep(-1, n), delta)), b = 0,
                 r = 0, l = matrix(rep(0, 3*n-1)), u = matrix(c(rep(meth_par[1], n-1), rep(meth_par[2], 2*n))), sigf = sigf,
                 maxiter = maxiter, margin = margin, bound = bound)
    #print("opt1")
    #print(head(opt1$xmin))
    #print(t(c(rep(0, n-1), rep(-1, n), delta))%*%opt1$xmin)
    #  print("opt2")
    #print(head(opt2@primal))
    # print(t(c(rep(0, n-1), rep(-1, n), delta)))
    # print(opt2@primal)
    # print(t(c(rep(0, n-1), rep(-1, n), delta))%*%opt2@primal)
    # if(abs(t(c(rep(0, n-1), rep(-1, n), delta))%*%opt2@primal) > 1) {
    #  browser()
    #}
    # print(H)
    #print(isSymmetric(H))
    # print(is.positive.semi.definite(round(H, sgf.sv)))
  }
  # Extract the parameters
  if (opt_alg == "quadprog") {
    alpha <- round(opt$xmin[1:(n-1)], sgf_sv)
    beta <- round(opt$xmin[n:(2*n - 1)], sgf_sv)
    beta.star <- round(opt$xmin[(2*n):(3*n - 1)], sgf_sv)
  } else {
    alpha <- round(opt@primal[1:(n-1)], sgf_sv)
    beta <- round(opt@primal[n:(2*n - 1)], sgf_sv)
    beta.star <- round(opt@primal[(2*n):(3*n - 1)], sgf_sv)
    #verification that the constraints are satisfied when 'ipop' is called
    if (abs(t(c(rep(0, n-1), rep(-1, n), delta))%*%opt@primal) > margin) {
      warning("Warning: constrains violated when calling 'ipop'.")
    }
  }
  #estimation of b_0
  beta.fact <- beta - delta*beta.star
  y_hat <- crossprod(matrix(alpha), crossprod(t(Dc), K)) + crossprod(matrix(beta.fact), K)
  b0 <- mean(Y - y_hat)
  hybo <- HybridObj(Alpha = alpha, Beta = beta, Betastar = beta.star, Delta = delta,
                    Xtrain = X, DifMat = dm, Kernel = Ker, OptMeth = opt_alg, b0 = b0)
  return(hybo)
}

#--- construct the hybridObj class
#' Constructs object of class \code{VB2FitObj}.\cr
#'
#' @title \code{HybridObj} (hybrid approach)
#' @param Alpha [\code{vector(1)}]\cr
#' A part of the solution of the quadratic optimization problem of interest.
#' @param Beta [\code{vector(1)}]\cr
#' A part of the solution of the quadratic optimization problem of interest.
#' @param Betastar [\code{vector(1)}]\cr
#' A part of the solution of the quadratic optimization problem of interest.
#' @param Delta [\code{vector(1)}]\cr
#' Vector of status 1 = no censored.
#' @param Xtrain [\code{matrix(1)}]\cr
#' Matrix of training data points.
#' @param DifMat [\code{\link{Diffmatrix}(1)}]\cr
#' Matrix used to maked differences between neighbor points.
#' @param Kernel [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#' @param OptMeth [\code{\link{Kernel}(1)}]\cr
#' Program used to solve the optimization problem.
#' @param b0 [\code{numeric(1)}]\cr
#' The estimated offset.
#'
#' @return [\code{HybridObj(1)}]
#' Object of class \code{Hybrid} containing elements:
#' \tabular{ll}{
#'    \code{Alpha} \tab Solution of the quadratic optimization problem, \cr
#'    \code{Xtrain} \tab Matrix of training points,\cr
#'    \code{DifMat} \tab Matrix used to made differences between neighbouring points. \cr
#'    \code{Kernel} \tab Kernel matrix, an object of class \code{Kernel},\cr
#'    \code{OptMeth} \tab Program used to solve the quadratic optimization problem.\cr
#'  }
#' @keywords internal
HybridObj <- function (Alpha = NULL, Beta = NULL, Betastar = NULL,
                      Delta = NULL, Xtrain = NULL, DifMat = NULL,
                      Kernel = NULL, OptMeth = NULL, b0 = NULL) {
  hybo <- VB1FitObj(Alpha = Alpha, Xtrain = Xtrain, DifMat = DifMat, Kernel = Kernel, OptMeth = OptMeth)
  hybo <- setBeta.HybridObj(hybo, beta = Beta)
  hybo <- setBetastar.HybridObj(hybo, betastar = Betastar)
  hybo <- setDelta.HybridObj(hybo, delta = Delta)
  hybo <- c(hybo, b0 = b0)
  class(hybo) <- append(class(hybo), "HybridObj")
  class(hybo) <- append(class(hybo), "RegFitObj")
  class(hybo) <- append(class(hybo), "VB1FitObj")
  return(hybo)
}

#' Creator of the generic mutator \code{setBetastar}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @param betastar [\code{vector}]\cr
#' New value.
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
setBetastar <- function (hybo, betastar) {
  UseMethod("setBetastar", hybo)
}
#' Creator of the generic mutator \code{setDelta}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @param delta [\code{vector(1)}]\cr
#' New value.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setDelta <- function(hybo, delta) {
  UseMethod("setDelta", hybo)
}
#' Default mutator of the field \code{Beta} of the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @param betastar [\code{vector(1)}]\cr
#' New value.
#'
#' @return [\code{Hybrid(1)}]
#' The object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setBetastar.default <- function(hybo, betastar) {
  return(hybo)
}
#' Default mutator of the field \code{delta} of the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @param delta [\code{vector(1)}]\cr
#' New value
#'
#' @return [\code{Hybrid(1)}]
#' The object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setDelta.default <- function(hybo, delta) {
  return(hybo)
}
#' Default mutator of the field \code{Beta} of the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param rfo [\code{Hybrid(1)}]\cr
#' Object of class \code{Hybrid} taken in the argument.
#' @param beta [\code{vector(1)}]\cr
#' Index of binary/categorial variables.
#'
#' @return [\code{Hybrid(1)}]
#' Modified version of the object taken in the argument.
#' @keywords internal
setBeta.HybridObj <- function(rfo, beta) {
  rfo$Beta <- beta
  return(rfo)
}
#' Default mutator of the field \code{Betastar} of the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object of class \code{Hybrid} taken in the argument.
#' @param betastar [\code{vector(1)}]\cr
#' New value.
#'
#' @return [\code{Hybrid(1)}]
#' Modified version of the object taken in the argument.
#' @keywords internal
setBetastar.HybridObj <- function(hybo, betastar) {
  hybo$Betarstar <- betastar
  return(hybo)
}
#' Creator of generic setor \code{setDelta}
#'
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object of class \code{Hybrid} taken in the argument.
#' @param delta [\code{vector(1)}]\cr
#' New value.
#' 
#' @return [\code{Hybrid(1)}]
#' Modified version of the object taken in the argument.
#' @keywords internal
setDelta.HybridObj <- function(hybo, delta) {
  hybo$Delta <- delta
  return(hybo)
}
#' Creator of the generic accessor \code{getBetastar}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBetastar <- function(hybo) {
  UseMethod("getBetastar", hybo)
}
#' Creator of the generic accessor \code{getDelta}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getDelta <- function(hybo) {
  UseMethod("getDelta", hybo)
}
#' Accessor for the field \code{Betastar} for the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBetastar.default <- function(hybo) {
  return(NULL)
}
#' Accessor for the field \code{Delta} for the object taken in an argument.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getDelta.default <- function(hybo) {
  return(NULL)
}
#' Creator of the generic accessor \code{Beta}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param rfo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @return [\code{vector(1)}]
#' Beta field of the object of class \code{Hybrid} taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBeta.HybridObj <- function(rfo) {
  return(rfo$Beta)
}
#' Creator of the generic accessor \code{Betastar}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in the argument.
#' @return [\code{vector(1)}]
#' Betastar field of the object of class \code{Hybrid} taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBetastar.HybridObj <- function(hybo) {
  return(hybo$Betarstar)
}
#' Creator of the generic accessor \code{Delta}.
#'
#'
#' @title \code{Hybrid} (hybrid approach)
#' @param hybo [\code{Hybrid(1)}]\cr
#' Object taken in an argument.
#' @return [\code{vector(1)}]
#' Delta field of the object of class \code{Hybrid} taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getDelta.HybridObj <- function(hybo) {
  return(hybo$Delta)
}
