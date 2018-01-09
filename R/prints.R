#' Prints  object of class \code{survivalsvm}.
#'
#'
#' @title print survivalsvm
#' @param x [\code{survivalsvm(1)}]\cr
#' Object of class \code{survivalsvm} to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#'
#' @export
#'
#' @author Cesaire J. K. Fouodo
print.survivalsvm <- function(x, ...) {
  # dispatching
  if (x$typeofsurvivalsvm == "regression") {
    result <- printRegFitObj(object = x)
  }
  if (x$typeofsurvivalsvm == "vanbelle1") {
    result <- printVB1FitObj(object = x)
  }
  if (x$typeofsurvivalsvm == "vanbelle2") {
    result <- printVB2FitObj(object = x)
  }
  if (x$typeofsurvivalsvm == "hybrid") {
    result <- printHybrid(object = x)
  }
}

#' Print of object of class \code{RegFitObj}. \code{RegFitObj} is the class of models fitted using the regression approach of survival support vector machines.
#'
#'
#' @title print survivalsvm
#' @param object [\code{survivalsvm(1)}]\cr
#' Object of class \code{survivalsvm} to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
printRegFitObj <- function(object, ...) {
  model <- object$model.fit
  cat("\n")
  cat("survivalsvm result\n\n")
  cat("Call:\n\n", deparse(object$call, width.cutoff = 500L), "\n\n")
  cat("Survival svm approach              :", object$typeofsurvivalsvm, "\n")
  cat("Type of Kernel                     :", getType(getKernel(model)), "\n")
  cat("Optimization solver used           :", getOptMeth(model), "\n")
  cat("Number of support vectors retained :", nrow(getSV.RegFitObj(model)), "\n")
  cat("survivalsvm version                :", 
      paste(object$package.version, collapse = "."), "\n")
}
#' Print of object of class \code{RegFitObj}. \code{VB1FitObj} is the class of models fitted using the ranking approach of survival support vector machines.
#'
#'
#' @title print survivalsvm
#' @param object [\code{survivalsvm(1)}]\cr
#' Object of class \code{survivalsvm} to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
printVB1FitObj <- function(object, ...) {
  model <- object$model.fit
  n.sv <- sum(getAlpha.VB1FitObj(model) != 0)
  cat("\n\n")
  cat("survivalsvm result\n\n")
  cat("Call:\n", deparse(object$call, width.cutoff = 500L), "\n\n")
  cat("Survival svm approach                  :", object$typeofsurvivalsvm, "\n")
  cat("Type of Kernel                         :", getType(getKernel(model)), "\n")
  cat("Method used to build 1NN differences   :", getType(getDifMat(model)), "\n")
  cat("Optimization solver used               :", getOptMeth(model), "\n")
  cat("Number of support vectors retained     :", n.sv, "\n")
  cat("survivalsvm version                    :", 
      paste(object$package.version, collapse = "."), "\n")
}


#' Print of object of class \code{RegFitObj}. \code{VB1FitObj} is the class of models fitted using the ranking approach of survival support vector machines.
#'
#'
#' @title print survivalsvm
#' @param object [\code{survivalsvm(1)}]\cr
#' Object \code{survivalsvm}  to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
printVB2FitObj <- function(object, ...) {
  model <- object$model.fit
  n.sv <- sum(getAlpha.VB1FitObj(model) != 0)
  cat("\n\n")
  cat("survivalsvm result\n\n")
  cat("Call:\n", deparse(object$call, width.cutoff = 500L), "\n\n")
  cat("Survival svm approach                  :", object$typeofsurvivalsvm, "\n")
  cat("Type of Kernel                         :", getType(getKernel(model)), "\n")
  cat("Method used to build 1NN differences   :", getType(getDifMat(model)), "\n")
  cat("Optimization solver used               :", getOptMeth(model), "\n")
  cat("Number of support vectors retained     :", n.sv, "\n")
  cat("survivalsvm version                    :", 
      paste(object$package.version, collapse = "."), "\n")
}


#' Print of object of class \code{Hybrid}. \code{Hybrid} is the class of models fitted using the hybrid approach of survival support vector machines.
#'
#'
#' @title print survivalsvm
#' @param object [\code{survivalsvm(1)}]\cr
#' Object \code{survivalsvm}  to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
printHybrid <- function(object, ...) {
  model <- object$model.fit
  n.sv <- sum(!(c(0, getAlpha.VB1FitObj(model)) != 0) | (getBeta.HybridObj(model) != 0) | (getBetastar.HybridObj(model) != 0))
  cat("\n\n")
  cat("survivalsvm result\n\n")
  cat("Call:\n", deparse(object$call, width.cutoff = 500L), "\n\n")
  cat("Survival svm approach                  :", object$typeofsurvivalsvm, "\n")
  cat("Type of kernel                         :", getType(getKernel(model)), "\n")
  cat("Method used to build 1NN differences   :", getType(getDifMat(model)), "\n")
  cat("Optimization solver used               :", getOptMeth(model), "\n")
  cat("Number of support vectors retained     :", n.sv, "\n")
  cat("survivalsvm version                    :", 
      paste(object$package.version, collapse = "."), "\n")
}

#' Print objects of class \code{survivalsvm}.
#'
#'
#' @title print survivalsvm
#' @param x [\code{survivalsvm(1)}]\cr
#' Object \code{survivalsvm}  to be printed.
#' @param ... [\code{any}]\cr
#' Further arguments passed to or from other methods.
#'
#' @author Cesaire J. K. Fouodo
#' @export
print.survivalsvmprediction <- function(x, ...) {
  cat("\n\n")
  cat("survivalsvm prediction\n\n")
  cat("Type of survivalsvm                      :", x$typeofsurvivalsvm, "\n")
  cat("Type of kernel                           :", x$typeofkernel, "\n")
  cat("Optimization solver used in model        :", x$opt.meth, "\n")
  if(x$typeofsurvivalsvm != "regression"){
  cat("Method used to build 1NN differences     :", x$diff.method,"\n")
  }
  cat("predicted risk ranks                     :", round(x$predicted[1:min(5, 
                                            length(x$predicted))], 2), "...\n")
  cat("survivalsvm version                      :", 
      paste(x$package.version, collapse = "."), "\n")
}
