#' Predictions of objects of class \code{survivalsvm}.
#'
#'
#' @title Suvirvalsvm predictions
#' @param object [\code{survivalsvm(1)}]\cr
#' Object of class \code{survivalsvm}, fitted with \code{\link{survivalsvm}}.
#' @param newdata [\code{data.frame(1)}]\cr
#' Data frame of observations.
#' @param subset [\code{vector(1)}]\cr
#' Indexes of data points of used to make the prediction.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#'
#' @return [\code{survivalsvmprediction(1)}]
#' Object of class \code{survivalsvmprediction}, with elements:
#' \tabular{ll}{
#'    \code{typeofsurvivalsvm} \tab Type of  \code{survivalsvm} object that is fitted in the model, \cr
#'    \code{typeofkernel} \tab type of kernel used to fit the model,\cr
#'    \code{parameterofkernel} \tab Kernel parameters used to fit the model,\cr
#'    \code{opt.meth} \tab solver used to fit the model,\cr
#'    \code{predicted} \tab values predicted.\cr
#'  }
#'
#' @author Cesaire J. K. Fouodo
#' @export
#'
#' @seealso \link{survivalsvm}
#' @examples
#' require(survival)
#' set.seed(123)
#' n <- nrow(veteran)
#' train.index <- sample(1:n, 0.7*n, replace = FALSE)
#' test.index <- setdiff(1:n, train.index)
#' survsvm.reg <- survivalsvm(Surv(veteran$diagtime, veteran$status) ~ .,
#'                            subset = train.index, data = veteran,
#'                            type = "regression", gamma.mu = 1,
#'                            opt.meth = "quadprog", kernel = "add_kernel")
#' pred.survsvm.reg <- predict(object = survsvm.reg, newdata = veteran, subset = test.index)
#' print(pred.survsvm.reg)
#' @importFrom utils packageVersion
predict.survivalsvm <- function(object, newdata, subset = NULL, ...) {
  if (!inherits(object, "survivalsvm")) {
    stop("Error: 'object' must be a survivalsvm object.")
  } else {
    if(is.null(object$var.names)) {
      stop("Error: Invalid 'survivalsvm' structure")
    }
    if (is.null(object$typeofsurvivalsvm)) {
      stop("Error: Invalid 'survivalsvm' structure")
    } else {
      if((object$typeofsurvivalsvm != "regression") & (object$typeofsurvivalsvm != "vanbelle1") &
         (object$typeofsurvivalsvm != "vanbelle2") & (object$typeofsurvivalsvm != "hybrid")) {
        stop("Error: Invalid 'survivalsvm' structure")
      }
    }
    if (is.null(getOptMeth(object$model.fit))) {
      stop("Error: Invalid 'survivalsvm' structure")
    } else {
      if(!(getOptMeth(object$model.fit) %in% c("quadprog", "ipop"))) {
        stop("Error: Invalid 'survivalsvm' structure")
      }
    }
    if (!inherits(object$model.fit, c("RegFitObj", "VB1FitObj", "VB2FitObj", "HybridObj"))) {
      stop("Error: 'survivalsvm' must inherit from 'RegFitObj', 'VB1FitObj', 'VB2FitObj' or 'HybridObj'.")
    }
  }
  kernel <- getKernel(object$model.fit)
  if (!inherits(kernel, "Kernel")) {
    stop("Error: Invalid Kernel.")
  } else {
    mat.kernel <- getMat(kernel)
    par.kernel <- getKernpar(kernel)
    bin.cat.kernel <- getBincat(kernel)
    type.kernel <- getType(kernel)
    if (is.null(mat.kernel) | is.null(par.kernel) |
       is.null(bin.cat.kernel) | is.null(type.kernel)) {
      stop("Error: Invalid Kernel.")
    }
  }
  if (inherits(object$model.fit, c("VB1FitObj", "VB2FitObj", "HybridObj"))) {
    mat.diff.obj <- getDifMat(object$model.fit)
    if (!inherits(mat.diff.obj, "Diffmatrix")) {
      stop("Error: Invalid matrix difference for neigborhood.")
    } else {
      mat.diff <- getMat(mat.diff.obj)
      mat.type <- getType(mat.diff.obj)
      if (is.null(mat.diff) | is.null(mat.type)) {
        stop("Error: Unknown Diffmatrix class.")
      }
    }
  }
  if (inherits(object$model.fit, "RegFitObj") & !inherits(object$model.fit, "HybridObj")) {
    beta <- getBeta(object$model.fit)
    sv <- getSV(object$model.fit)
    b0 <- getb0.RegFitObj(object$model.fit)
    if (is.null(beta) | is.null(sv) | is.null(b0)) {
      stop("Error: Invalid 'RegFitObj' object.")
    }
  }
  if (inherits(object$model.fit, c("VB1FitObj", "VB2FitObj"))) {
    # Alpha = NULL, Xtrain = NULL, DifMat = NULL, Kernel = NULL
    alpha <- getAlpha(object$model.fit)
    train <- getXtrain(object$model.fit)
    if(is.null(alpha) | is.null(train)) {
      stop("Error: Unknown structure.")
    }
  }
  if (inherits(object$model.fit, "HybridObj")){
    alpha <- getAlpha(object$model.fit)
    beta <- getBeta(object$model.fit)
    betastar <- getBetastar(object$model.fit)
    delta <- getDelta(object$model.fit)
    train <- getXtrain(object$model.fit)
    b0 <- getb0(object$model.fit)
    if (is.null(alpha) | is.null(beta) |
       is.null(betastar) | is.null(delta) |
       is.null(train) | is.null(b0)) {
      stop("Error: Invalid 'HybridObj' object.")
    }
  }
  # Test on newdata
  if (!inherits(newdata, "data.frame")) {
    stop("Error: 'newdata' must be an object from the class 'data.frame'.")
  }
  if (!all(object$var.names %in% names(newdata))) {
    stop("Error: names of independing variables in 'newdata'  do not match with those in the 'object'.")
  }
  X_pred <- if(is.null(subset)){
    data.matrix(newdata[, object$var.names])
  } else {
    data.matrix(newdata[subset, object$var.names])
  }
  X_pred <- X_pred[stats::complete.cases(X_pred),]
  if (ncol(X_pred) == 0) {
    stop("No observation in the data frame given in argument.")
  }
  # dispatching
  if (object$typeofsurvivalsvm == "regression") {
    result <- predictRegFitObj(object = object$model.fit, X_pred = X_pred)
  }
  if (object$typeofsurvivalsvm == "vanbelle1") {
    result <- predictVB1FitObj(object = object$model.fit, X_pred = X_pred)
  }
  if (object$typeofsurvivalsvm == "vanbelle2") {
    result <- predictVB2FitObj(object = object$model.fit, X_pred = X_pred)
  }
  if (object$typeofsurvivalsvm == "hybrid") {
    result <- predictHybrid(object = object$model.fit, X_pred = X_pred)
  }
  result$package.version <- unlist(packageVersion("survivalsvm"))
  return(result)
}

#' Predictions based on model fitted using the regression approach of survival support vector machines.
#'
#'
#' @title Survivalsvm predictions
#' @param object [\code{\link{survivalsvm}(1)}]\cr
#' makes predictions using a survivalsvm object fitted with \link{survivalsvm}.
#' @param X_pred [\code{matrix(1)}]\cr
#' matrix of data points of interest.
#'
#' @return [\code{survivalsvmprediction}(1)]
#' object of class survivalsvmprediction, with elements:
#' \tabular{ll}{
#'    \code{typeofsurvivalsvm} \tab type of the survivalsvm object that is fitted in model, \cr
#'    \code{typeofkernel} \tab type of kernel used to fit the model,\cr
#'    \code{parameterofkernel} \tab parameters of kernel that used to fit the model,\cr
#'    \code{opt.meth} \tab program used to fit the model,\cr
#'    \code{predicted} \tab values predicted.\cr
#'  }
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
predictRegFitObj <- function(object, X_pred){
  Xtrain <- getSV.RegFitObj(object)
  ker <- getKernel.RegFitObj(object)
  kernel_type <- getType.Kernel(ker)
  kernel_pars <- getKernpar.Kernel(ker)
  bin_cat <- getBincat.Kernel(ker)
  K.pred <- kernelMatrix(Xtrain = getSV.RegFitObj(object),
                         kernel_type = kernel_type,
                         kernel_pars = kernel_pars,
                         Xt = X_pred,
                         bin_cat = bin_cat)
  beta <- getBeta.RegFitObj(object)
  b0 <- getb0.RegFitObj(object)
  y_hat <- crossprod(matrix(beta), getMat.Kernel(K.pred)) + b0
  result <- list(
    typeofsurvivalsvm = "regression",
    typeofkernel = kernel_type,
    parameterofkernel = kernel_pars,
    opt.meth = getOptMeth(object),
    predicted = y_hat
  )
  class(result) <- "survivalsvmprediction"
  return(result)
}

#' Predictions based on model fitted using the ranking approach of survival support vector machines.
#'
#'
#' @title Survivalsvm predictions
#' @param object [\code{\link{survivalsvm}(1)}]\cr
#' survivalsvm object, fitted with \link{survivalsvm}.
#' @param X_pred [\code{matrix}]\cr
#' matrix of data points of interest.
#'
#' @return object of class survivalsvmprediction, with elements:
#' \tabular{ll}{
#'    \code{typeofsurvivalsvm} \tab type of the survivalsvm object that is fitted in model, \cr
#'    \code{typeofkernel} \tab type of kernel used to fit the model,\cr
#'    \code{parameterofkernel} \tab parameters of kernel that used to fit the model,\cr
#'    \code{opt.meth} \tab program used to fit the model,\cr
#'    \code{predicted} \tab values predicted.\cr
#'  }
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
predictVB1FitObj <- function(object, X_pred){
  train <- getXtrain.VB1FitObj(object)
  kerntype <- getType.Kernel(getKernel.VB1FitObj(object))
  kernpar <- getKernpar.Kernel(getKernel.VB1FitObj(object))
  bin_cat <- getBincat.Kernel(getKernel.VB1FitObj(object))
  Dc <- getMat.Diffmatrix(getDifMat.VB1FitObj(object))
  alpha <- getAlpha.VB1FitObj(object)
  Kern.pred <- kernelMatrix(Xtrain = train,
                            kernel_type = kerntype,
                            kernel_pars = kernpar,
                            Xt = X_pred, bin_cat = bin_cat)
  K.pred <- getMat.Kernel(Kern.pred)
  y_hat <- crossprod(matrix(alpha), crossprod(t(Dc), K.pred))
  result <- list(
    typeofsurvivalsvm = "vanbelle1",
    typeofkernel = kerntype,
    parameterofkernel = kernpar,
    opt.meth = getOptMeth(object),
    diff.method = getType.Diffmatrix(getDifMat.VB1FitObj(object)),
    predicted = y_hat
  )
  class(result) <- "survivalsvmprediction"
  return(result)
}

#' Predictions based on model fitted using the ranking approach of survival support vector machines.
#'
#'
#' @title Survivalsvm predictions
#' @param object [\code{\link{survivalsvm}(1)}]\cr
#' Object of class \code{survivalsvm}, fitted with \link{survivalsvm}.
#' @param X_pred [\code{matrix}]\cr
#' Matrix of data points used to make the prediction.
#'
#' @return object of class survivalsvmprediction, with elements:
#' \tabular{ll}{
#'    \code{typeofsurvivalsvm} \tab Type of the survivalsvm object that is fitted in model, \cr
#'    \code{typeofkernel} \tab type of kernel used to fit the model,\cr
#'    \code{parameterofkernel} \tab parameters of kernel that used to fit the model,\cr
#'    \code{opt.meth} \tab program used to fit the model,\cr
#'    \code{predicted} \tab values predicted.\cr
#'  }
#' @keywords internal
predictVB2FitObj <- function(object, X_pred){
  train <- getXtrain.VB1FitObj(object)
  kerntype <- getType.Kernel(getKernel.VB1FitObj(object))
  kernpar <- getKernpar.Kernel(getKernel.VB1FitObj(object))
  bin_cat <- getBincat.Kernel(getKernel.VB1FitObj(object))
  Dc <- getMat.Diffmatrix(getDifMat.VB1FitObj(object))
  alpha <- getAlpha.VB1FitObj(object)
  Kern.pred <- kernelMatrix(Xtrain = train,
                            kernel_type = kerntype,
                            kernel_pars = kernpar,
                            Xt = X_pred, bin_cat = bin_cat)
  K.pred <- getMat.Kernel(Kern.pred)
  y_hat <- crossprod(matrix(alpha), crossprod(t(Dc), K.pred))
  result <- list(
    typeofsurvivalsvm = "vanbelle2",
    typeofkernel = kerntype,
    parameterofkernel = kernpar,
    opt.meth = getOptMeth(object),
    diff.method = getType.Diffmatrix(getDifMat.VB1FitObj(object)),
    predicted = y_hat
  )
  class(result) <- "survivalsvmprediction"
  return(result)
}

#' Maker of predictions based on model fitted using the hybrid approach of survival support vector machines.
#'
#'
#' @title Survivalsvm predictions
#' @param object [\code{\link{survivalsvm}(1)}]\cr
#' survivalsvm object, fitted with \link{survivalsvm}.
#' @param X_pred [\code{matrix}]\cr
#' matrix of data points of interest.
#'
#' @return object of class survivalsvmprediction, with elements:
#' \tabular{ll}{
#'    \code{typeofsurvivalsvm} \tab type of the survivalsvm object that is fitted in model, \cr
#'    \code{typeofkernel} \tab type of kernel used to fit the model,\cr
#'    \code{parameterofkernel} \tab parameters of kernel used to fit the model,\cr
#'    \code{opt.meth} \tab program used to fit the model,\cr
#'    \code{predicted} \tab values predicted.\cr
#'  }
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
predictHybrid <- function(object, X_pred){
  train <- getXtrain.VB1FitObj(object)
  K <- getKernel.VB1FitObj(object)
  kernel_type <- getType.Kernel(K)
  kernel_pars <- getKernpar.Kernel(K)
  bin_cat <- getBincat.Kernel(K)
  alpha <- getAlpha.VB1FitObj(object)
  beta <- getBeta.HybridObj(object)
  beta.star <- getBetastar.HybridObj(object)
  delta <- getDelta.HybridObj(object)
  beta.fact <- beta - delta*beta.star
  Dc <- getMat.Diffmatrix(getDifMat.VB1FitObj(object))
  b0 <- getb0.RegFitObj(object)
  Kern <- kernelMatrix(Xtrain = train,
                       kernel_type = kernel_type,
                       kernel_pars = kernel_pars,
                       Xt = X_pred, bin_cat = bin_cat)
  Kern.pred <- getMat.Kernel(Kern)
  y_hat <- crossprod(matrix(alpha), crossprod(t(Dc), Kern.pred)) + crossprod(matrix(beta.fact), Kern.pred) + b0
  result <- list(
    typeofsurvivalsvm = "hybrid",
    typeofkernel = kernel_type,
    parameterofkernel = kernel_pars,
    opt.meth = getOptMeth(object),
    diff.method = getType.Diffmatrix(getDifMat.VB1FitObj(object)),
    predicted = y_hat
  )
  class(result) <- "survivalsvmprediction"
  return(result)
}
