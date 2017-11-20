# function regFit fit the regression approach of svm for survival time analysis
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
#   kernel_matrix for more details.
# @param quadprog [character(1)]
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
# @return [RegFitObj(1)]
#              beta.fact [vector(1)] estimated factors
#              sv [interger(1)] number of support vectors
#              kernel_type [character] kernel used during fitting phase
#              kernel_pars [vector(1)] parameters used to construict the kernel matrix
#              b0 [interger(1)] the estimated bias
#---------------------------------------------------------------------------------------------
#' The function \code{regFit} fits a \code{survivalsvm} model based on the regression approach.
#'
#'
#' @title survivalsvm (regression approach)
#' @param X [\code{matrix(1)}]\cr
#' design matrix.
#' @param Y [\code{vector(1)}]\cr
#' vector of survival times.
#' @param delta [\code{vector(1)}]\cr
#' vector of status: 0 if censored and 1 else.
#' @param meth_par [\code{numeric(1)}]\cr
#' parameter of regularization.
#' @param kernel_type [\code{character(1)}]\cr
#' type of the kernel.
#' @param kernel_pars [\code{vector(1)}]\cr
#' parameter of kernel.
#' @param bin_cat [\code{vector(1)}]\cr
#' indexes of binary/categorial variables.
#' @param opt_alg [\code{character(1)}]\cr
#' program used to solve the optimization problem. This most be one of the two possibilities \code{\link{quadprog}} or \code{\link{ipop}}.
#' @param sgf_sv [\code{integer(1)}]\cr
#' number of digits to be retained in the solution.
#' @param sigf [\code{integer(1)}]\cr
#' used by \code{ipop}. See \code{\link{ipop}} for more details.
#' @param maxiter [\code{integer(1)}]\cr
#' used by \code{ipop}. See \code{\link{ipop}} for more details.
#' @param margin [\code{numeric(1)}]\cr
#' used by \code{ipop}. See \code{\link{ipop}} for more details.
#' @param bound [\code{numeric(1)}]\cr
#' used by \code{ipop}. See \code{\link{ipop}} for more details.
#' @param eig.tol [\code{numeric(1)}]\cr
#' used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#' @param conv.tol [\code{numeric(1)}]\cr
#' used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#' @param posd.tol [\code{numeric(1)}]\cr
#' used by \code{nearPD} for adjusting positive definiteness. See \code{\link{nearPD}} for detail.
#'
#' @export
#' @return [\code{\link{RegFitObj}(1)}]
#' object of class \code{RegFitObj} containing elements:
#' \tabular{ll}{
#'    \code{Beta} \tab solution of the quadratic optimization problem, \cr
#'    \code{SV} \tab support vector machines,\cr
#'    \code{Kernel} \tab kernel matrix, an object of class \code{Kernel},\cr
#'    \code{b0} \tab estimated offset,\cr
#'    \code{OptMeth} \tab program used to solve the quadratic optimization problem.\cr
#'  }
#'
#' @author Cesaire J. K. Fouodo
#' @keywords internal
regFit <- function(X, Y, delta, meth_par = 1, kernel_type = "lin_kernel",
                   kernel_pars = NA, bin_cat = integer(0), opt_alg = "quadprog",
                   sgf_sv = 5, sigf = 7, maxiter = 20, margin = 0.05, bound = 10,
                   eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08){
  # Implementation of the regression approach
  # Use censored failure time data
  # event(not censored): delta=1;  censored: delta=0
  if (!(opt_alg %in% c("quadprog", "ipop"))) {
    stop("'opt_alg' must be either 'quadprog' or 'ipop'")
  }
  n <- nrow(X)
  d <- ncol(X)
  if (is.na(kernel_pars) & !(kernel_type == "add_kernel" || kernel_type == "lin_kernel")) {
    kernel_pars <- rep(1/ncol(X), ncol(X))
  }
  # Build the Kernel matrix
  ker <- kernelMatrix(Xtrain = X, kernel_type = kernel_type, kernel_pars = kernel_pars, bin_cat = bin_cat)
  K <- getMat.Kernel(ker)
  # Solve dual problem by quadprog
  D1 <- cbind(K, -delta * K)
  D2 <- cbind(-delta * K, K)
  D <- rbind(D1, D2)
  opt <- if (opt_alg == "quadprog") {
    quadprog::quadprog(C =  as.matrix(Matrix::nearPD(D, eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)$mat),
                     d = c(-Y, delta*Y), Aeq = c(rep(-1, n), delta),
             beq = 0, A = -diag(2*n), b = rep(0, 2*n), lb = 0, ub = meth_par)
  } else{
    kernlab::ipop(H = D, c = c(-Y, delta*Y), A = t(c(rep(-1, n), delta)), b = 0, l = matrix(0, 2*n),
         u = matrix(meth_par, 2*n), r = 0, sigf = sigf, maxiter = maxiter, margin = margin, bound = bound)
  }
  #estimation of b_0
  betapar <- if (opt_alg == "quadprog"){
    cbind(round(opt$xmin[1:n], sgf_sv), round(opt$xmin[(n + 1):(2 * n)], sgf_sv))
  } else {
    cbind(round(opt@primal[1:n], sgf_sv), round(opt@primal[(n + 1):(2 * n)], sgf_sv))
  }
  i.sv <- which(!(rowSums(betapar) == 0))
  if(length(i.sv) == 0) {
    stop("Inconsistent solutions when trying optimization with the given parameters.")
  }
  betapar <- matrix(betapar[i.sv,], ncol = 2, byrow = FALSE)
  sv <- matrix(X[i.sv, ], byrow = FALSE, ncol = ncol(X))
  Ker.sv <- kernelMatrix(Xtrain = sv, kernel_type = kernel_type, kernel_pars = kernel_pars, bin_cat = bin_cat)
  K.sv <- getMat.Kernel(Ker.sv)
  beta.fact <- betapar[,1] - delta[i.sv] * betapar[,2]
  y_hat <- crossprod(matrix(beta.fact), K.sv)
  b0 <- mean(Y[i.sv] - y_hat)
  rfo <- RegFitObj(Beta = beta.fact, SV = sv, Kernel = ker, b0 = b0, OptMeth = opt_alg)
  return(rfo)
}

# function regPredict uses a model fitted with the regression approach to perform
# @param X_pred [matrix(1)]
#   matrix of observations
# @param modelpar[list(1)]
#   list of returned by the regFit function. See regFit function for more details
# @return [vector(1)]
#   Vector of predictions using the regression approach
#--------------------------------------------------------------------------------
## Makes predictions using using a model fitted with the regression survival support vector machines approach.
##
## @param X_pred data ponts of interest.
## @param model object of class \code{regFit}.
##
## @return vector of predicted values.
## @export
# regPredict <- function(X_pred, model){
#   Xtrain <- getSV.RegFitObj(model)
#   ker <- getKernel.RegFitObj(model)
#   kernel_type <- getType.Kernel(ker)
#   kernel_pars <- getKernpar.Kernel(ker)
#   bin_cat <- getBincat.Kernel(ker)
#   K.pred <- kernelMatrix(Xtrain = getSV.RegFitObj(model),
#                           kernel_type = kernel_type,
#                           kernel_pars = kernel_pars,
#                           Xt = X_pred,
#                           bin_cat = bin_cat)
#   beta <- getBeta.RegFitObj(model)
#   b0 <- getb0.RegFitObj(model)
#   y_hat <- crossprod(matrix(beta), getMat.Kernel(K.pred)) + b0
#   return(y_hat)
# }

# >>>> contuction of RegFitObj class
  #--- the constructor
#' Constructs object of class \code{RegFitObj}.
#'
#' @title survivalsvm (regression approach)
#' @param Beta [\code{vector(1)}]\cr
#' solution of the quadratic optimization problem of interest
#' @param SV [\code{matrix(1)}]\cr
#' support vector machines.
#' @param Kernel [\code{\link{Kernel}(1)}]\cr
#' object of class \code{Kernel}.
#' @param OptMeth [\code{character(1)}]\cr
#' program used to solve the optimization problem.
#' @param b0 [\code{numeric(1)}]\cr
#' the estimated offset.
#'
#' @return [\code{RegFitObj}(1)]
#' object of class \code{RegFitObj} containing elements:
#'
#' @keywords internal
  RegFitObj <- function(Beta = NULL, SV = NULL, Kernel = NULL, OptMeth = NULL, b0 = NULL) {
    rfo <- list(Beta = Beta, SV = SV, Kernel = Kernel, OptMeth = OptMeth, b0 = b0)
    class(rfo) <- "RegFitObj"
    return(rfo)
  }
#--- Generic mutator method for the Diffmatrix type
#--------------------------------------------------
#' Creator of the generic mutator \code{setBeta}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param beta [\code{vector(1)}]\cr
#' new value.
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
  setBeta <- function(rfo, beta) {
    UseMethod("setBeta", rfo)
  }
#' Creator of the generic mutator \code{setSV}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param sv [\code{matrix(1)}]\cr
#' new value.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setSV <- function(rfo, sv) {
  UseMethod("setSV", rfo)
}
#' Creator of the generic mutator \code{setKernel}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param kernel [\code{\link{Kernel}(1)}]\cr
#' new value.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setKernel <- function(rfo, kernel) {
  UseMethod("setKernel", rfo)
}
#' Creator of the generic mutator \code{setb0}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param b0 [\code{numeric(1)}]\cr
#' new value.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setb0 <- function(rfo, b0) {
  UseMethod("setb0", rfo)
}
#' Creator of the generic mutator \code{setOptMeth}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param optmeth [\code{character(1)}]\cr
#' new value.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setOptMeth <- function(rfo, optmeth) {
  UseMethod("setOptMeth", rfo)
}
#' Creator of the generic accessor \code{getBeta}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBeta <- function(rfo) {
  UseMethod("getBeta", rfo)
}
#' Creator of the generic accessor \code{getSV}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getSV <- function(rfo) {
  UseMethod("getSV", rfo)
}
#' Creator of the generic accessor \code{getKernel}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernel <- function(rfo) {
  UseMethod("getKernel", rfo)
}
#' Creator of the generic accessor \code{getOptMeth}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getOptMeth <- function(rfo) {
  UseMethod("getOptMeth", rfo)
}
#' Creator of the generic accessor \code{getb0}.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getb0 <- function(rfo) {
  UseMethod("getb0", rfo)
}
#>>>> Mutators for the kernel object
#' Default mutator of the field \code{Beta} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param beta [\code{vector(1)}]\cr
#' new offset.
#'
#' @return the object taken in an argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setBeta.default <- function(rfo, beta) {
  return(rfo)
}
#' Default mutator of the field \code{SV} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param sv [\code{matrix(1)}]\cr
#' new support vectors.
#'
#' @return [\code{RegFitObj}(1)]
#' the object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setSV.default <- function(rfo, sv) {
  return(rfo)
}

#' setKernel.default
#'
#'
#' @title \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param Kernel [\code{\link{Kernel}(1)}]\cr
#' new object of class \code{Kernel}.
#'
#' @return [\code{RegFitObj}(1)]
#' modified object.
#' @keywords internal
setKernel.default <- function(rfo, Kernel) {
  return(rfo)
}
#' Default mutator of the field \code{OptMeth} of the object taken in an argument.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param optmeth [\code{character(1)}]\cr
#' new name the of optimization program.
#'
#' @return [\code{RegFitObj}(1)]
#' the object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setOptMeth.default <- function(rfo, optmeth) {
  return(rfo)
}
#' Default mutator of the field \code{b0} of the object taken in an argument.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @param b0 [\code{numeric(1)}]\cr
#' new offset value.
#'
#' @return [\code{RegFitObj}(1)]
#' the object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setb0.default <- function(rfo, b0) {
  return(rfo)
}
#' Accessor for the field \code{Beta} for the object taken in an argument.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBeta.default <- function(rfo) {
  return(NULL)
}
#' Accessor for the field \code{SV} for the object taken in an argument.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getSV.default <- function(rfo) {
  return(NULL)
}

#' Accessor for the field \code{Kernel} for the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernel.default <- function(rfo) {
  return(NULL)
}
#' Accessor for the field \code{OptMeth} for the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getOptMeth.default <- function(rfo) {
  return(NULL)
}
#--- specific mutators
#' Default mutator of the field \code{Beta} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object of class \code{RegFitObj} taken in the argument.
#' @param beta [\code{vector(1)}]\cr
#' solutions of quadratic optimization problem.
#'
#' @return [\code{RegFitObj}(1)]
#' modified version of the object taken in the argument.
#' @keywords internal
setBeta.RegFitObj <- function(rfo, beta) {
  rfo$Beta <- beta
  return(rfo)
}
#--- specific mutators
#' Default mutator of the field \code{SV} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object of class \code{RegFitObj} taken in the argument.
#' @param sv [\code{matrix(1)}]\cr
#' support vectors.
#'
#' @return [\code{RegFitObj}(1)]
#' modified version of the object taken in the argument.
#' @keywords internal
setSV.RegFitObj <- function(rfo, sv) {
  rfo$SV <- sv
  return(rfo)
}
#--- specific mutators
#' Default mutator of the field \code{Kernel} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object of class \code{RegFitObj} taken in the argument.
#' @param kernel [\code{\link{Kernel}(1)}]\cr
#' index of binary/categorial variables.
#'
#' @return [\code{RegFitObj}(1)]
#' modified version of the object taken in the argument.
#' @keywords internal
setKernel.RegFitObj <- function(rfo, kernel) {
  rfo$Kernel <- kernel
  return(rfo)
}
#--- specific mutators
#' Default mutator of the field \code{OptMeth} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object of class \code{RegFitObj} taken in the argument.
#' @param optmeth [\code{character(1)}]\cr
#' names the solver.
#'
#' @return [\code{RegFitObj}(1)]
#' modified version of the object taken in the argument.
#' @keywords internal
setOptMeth.RegFitObj <- function(rfo, optmeth) {
  rfo$OptMeth <- optmeth
  return(rfo)
}
#--- specific mutators
#' Default mutator of the field \code{b0} of the object taken in an argument.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object of class \code{RegFitObj} taken in the argument.
#' @param b0 [\code{numeric(1)}]\cr
#' new offset.
#'
#' @return [\code{RegFitObj}(1)]
#' modified version of the object taken in the argument.
#' @keywords internal
setb0.RegFitObj <- function(rfo, b0) {
  rfo$b0 <- b0
  return(rfo)
}
#' Creator of the generic accessor \code{Beta}.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @return [\code{vector(1)}]
#' Beta field of the object of class \code{RegFitObj} taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBeta.RegFitObj <- function(rfo) {
  return(rfo$Beta)
}
#' Creator of the generic accessor \code{SV}.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @return [\code{matrix}]
#' the matrix of support vectors.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getSV.RegFitObj <- function(rfo) {
  return(rfo$SV)
}
#' Creator of the generic accessor \code{Kernel}.
#'
#'
#' @title \code{RegFitObj}
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @return [\code{Kernel}]
#' kernel.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernel.RegFitObj <- function(rfo) {
  return(rfo$Kernel)
}
#' Creator of the generic accessor \code{OptMeth}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @return [character(1)]
#' the named of the optimization program.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getOptMeth.RegFitObj <- function(rfo) {
  return(rfo$OptMeth)
}
#' Creator of the generic accessor \code{b0}.
#'
#'
#' @title Class \code{RegFitObj} (regression approach)
#' @param rfo [\code{RegFitObj}(1)]\cr
#' object taken in the argument.
#' @return [\code{numeric(1)}]
#' the offset
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getb0.RegFitObj <- function(rfo) {
  return(rfo$b0)
}
