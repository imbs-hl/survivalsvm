# function kernelMatrix to perform the kernel operation implicitly, without using of the feature transformation
#
# @param Xtrain [matrix(1)]
#   Data set for fitting.
# @param kernel_type [character(1)]
#   Name of the kernel typ to be used
# @param Xt [matrix(1)]
#   Matrix on which the mapping will be done. It is use to be the description matrix of
#   indidiums, for which the failure times are needed to be predict
# @param bin_cat [vector(1)]
#   index of binar or categorial columns in the input matrix
# @return [matrix(1)].
# The symmetric Kernel matrix
#----------------------------------------------------------------------------------------
#' computes the kernel matrix for elements taken in an argument.
#'
#' \code{Kernel}
#'
#'
#' @title Kernel
#' @param Xtrain [\code{matrix(1)}]\cr
#' Matrix of training data points.
#' @param kernel_type [\code{character(1)}]\cr
#' Type of kernel that is required.
#' @param kernel_pars [\code{vector(1)}]\cr
#' Parameters of kernels.
#' @param Xt [\code{matrix(1)}]\cr
#' Matrix of data points to be mapped.
#' @param bin_cat indexes of binary/categorical variables.
#'
#' @return [\code{\link{Kernel}(1)}]
#' Object of class Kernel, with elements:
#' \tabular{ll}{
#'    \code{Type} \tab type of kernel, \cr
#'    \code{Mat} \tab matrix used to compute differences between comparable data points. \cr
#'  }
#'
#' @author Cesaire J. K. Fouodo
kernelMatrix <- function (Xtrain, kernel_type = "lin_kernel", kernel_pars, Xt = NULL, bin_cat = integer(0)) {
  nb_data <- nrow(Xtrain)
  d <- ncol(Xtrain)
  if (tolower(kernel_type) == "rbf_kernel") {
    if (is.null(Xt)) {
      XXh <- crossprod(t(rowSums(Xtrain^2)), t(rep(1, nb_data)))
      omega <- XXh + t(XXh) - (2 * tcrossprod(Xtrain, Xtrain))
      omega <- exp(-omega / (2 * kernel_pars[1]))
    } else {
      XXh1 <- crossprod(t(rowSums(Xtrain^2)),t(rep(1, nrow(Xt))))
      XXh2 <- crossprod(t(rowSums(Xt^2)), t(rep(1, nb_data)))
      omega <- XXh1 + t(XXh2) - 2 * tcrossprod(Xtrain, Xt)
      omega <- exp(-omega / (2*kernel_pars[1]))
    }
  } else {
    if(tolower(kernel_type) == "rbf4_kernel") {
      if (is.null(Xt)) {
        XXh <- crossprod(t(rowSums(Xtrain^2)), t(rep(1, nb_data)))
        omega <- XXh + t(XXh) - (2 * tcrossprod(Xtrain, Xtrain))
        omega <- 0.5*(3 - t(apply(omega, 2, function(i) i/ kernel_pars))) * as.vector(exp(-omega / (2 * kernel_pars[1])))
      } else {
        XXh1 <- crossprod(t(rowSums(Xtrain^2)),  t(rep(1, nrow(Xt))))
        XXh2 <- crossprod(t(rowSums(Xt^2)), t(rep(1, nb_data)))
        omega <- XXh1 + t(XXh2) - 2 * tcrossprod(Xtrain, Xt)
        omega <- 0.5*(3 - omega / kernel_pars) * as.vector(exp(-omega / (2 * kernel_pars[1])))
      }
    } else {
      if (tolower(kernel_type) == "lin_kernel") {
        if (is.null(Xt)) {
          omega <- tcrossprod(Xtrain, Xtrain)
        } else {
          omega <- tcrossprod(Xtrain, Xt)
        }
      } else {
        if (tolower(kernel_type) == "poly_kernel") {
          if(is.null(Xt)){
            omega <- (crossprod(Xtrain, Xtrain) + kernel_pars[1])^kernel_pars[2]
          } else {
            omega <- (crossprod(Xtrain, Xt) + kernel_pars[1])^kernel_pars[2]
          }
        } else {
          if (tolower(kernel_type) == "add_kernel") {
            if(is.null(Xt)){
              Xt <- Xtrain
            }
            d <- ncol(Xtrain)
            num.ord <- setdiff(1:d, bin_cat)
            X.min <- apply(Xtrain[,num.ord], 2, min)
            X.max <- apply(Xtrain[,num.ord], 2, max)
            c.diff <- X.max - X.min
            if(any(c.diff == 0)){
              stop("additiv kernel can not be applied on constant column")
            }
            nr <- nrow(Xtrain)
            sub.X <- sapply(1:nr, function(i) {
              kp <- matrix(prod(dim(Xt)), ncol = ncol(Xt), nrow = nrow(Xt))
              kp[,num.ord] <- t((c.diff - abs(t(Xt[,num.ord]) - Xtrain[i, num.ord])) / c.diff)
              if (!(length(bin_cat) == 0)) {
                kp[, bin_cat] <- t(t(Xt[, bin_cat]) == Xtrain[i, bin_cat])
              }
              rs <- rowSums(kp)
            }, simplify = FALSE)
            omega <- do.call("rbind", sub.X)
          } else {
            stop("Unkonwn Kernel")
          }
        }
      }
    }
  }
  kern <- Kernel(Type = kernel_type, Mat = omega, Kernpar = kernel_pars, bincat = bin_cat)
  return(kern)
}

#--- Construct the kernel classes
#' Constructor of objects of class \code{Kernel}.
#'
#'
#' @title \code{Kernel}
#' @param Type  [\code{character}]\cr
#' Type of kernel. Must be chosen from the following strings: \code{"lin_kernel"}, \code{"add_kernel"},
#' \code{"rbf_kernel"}, \code{"rbf4_kernel"} or \code{"poly_kernel"}.
#' @param Mat [\code{matrix(1)}]\cr
#' Kernel matrix.
#' @param Kernpar [\code{vector(1)}]\cr
#' Kernel parameters, when required.
#' @param bincat [\code{vector(1)}]\cr
#' Index of binary/categorical variables.
#'
#' @return [\code{\link{Kernel}(1)}]
#' Object of class Kernel, with elements:
#' \tabular{ll}{
#'    \code{Type} \tab type of kernel, \cr
#'    \code{Mat} \tab matrix used to perform differences between comparable data points. \cr
#'  }
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
Kernel <- function(Type = NULL, Mat = NULL, Kernpar = NULL, bincat = NULL){
  k <- list(Type = Type, Mat = Mat, Kernpar = Kernpar, bincat = bincat)
  class(k) <- append(class(k), "Kernel")
  return(k)
}

#--- Generic mutator method for the kernel type
#----------------------------------------------
#' Creator of the generic mutator \code{setType}.
#'
#'
#' @title Class \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param value [\code{character(1)}]\cr
#' New type.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setType <- function(kern, value) {
  UseMethod("setType", kern)
}

#>>>> Mutators for the kernel object

#--- Default mutator method for the kernel type
#' Mutator of the field \code{Type} of the object taken in an argument.
#'
#'
#' @title Class \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param value [\code{character(1)}]\cr
#' New type.
#'
#' @return Object taken in the argument.
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setType.default <- function(kern, value) {
  return(kern)
}

#--- Specific mutator method for the kernel type
#' Mutator of the field \code{Type} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param kerntype [\code{character(1)}]\cr
#' Kerntype new type.
#'
#' @return Object of class \code{Kernel} with elements:
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setType.Kernel <- function(kern, kerntype) {
  kern$Type <- kerntype
  return(kern)
}

#--- Generic mutator method for the kernel matrix
#' Creator of generic mutator \code{setMat}.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param mat [\code{matrix(1)}]\cr
#' New matrix.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setMat <- function(kern, mat) {
  UseMethod("setMat", kern)
}
#--- Default mutator method for the kernel matrix
#' Mutator of the field \code{Mat} of the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param kernmatrix [\code{matrix(1)}]\cr
#' New matrix.
#'
#' @return The object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setMatrix.default <- function(kern, kernmatrix) {
  return(kern)
}

#--- Specific mutator method for the kernel matrix
#' Mutator of the field \code{Mat} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel} taken in the argument.
#' @param kernmat [\code{matrix(1)}]\cr
#' New kernel matrix.
#'
#' @return [\code{\link{Kernel}(1)}]
#' Object of class \code{Kernel} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab type of kernel, \cr
#'    \code{Mat} \tab kernel matrix, \cr
#'    \code{Kernpar} \tab parameters of kernel, when required, \cr
#'    \code{bincat} \tab index of binary/categorical variables, when required.
#'  }
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setMat.Kernel<-function(kern, kernmat) {
  kern$Mat <- kernmat
  return(kern)
}

#--- Generic mutator method for the kernel parameter
#---------------------------------------------------
#' Default mutator of the field \code{Kernpar} of the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel} taken in the argument.
#' @param kernpar [\code{vector(1)}]\cr
#' New kernel parameters.
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setKernpar <- function(kern, kernpar) {
  UseMethod("setMat", kern)
}
#--- Default mutator method for the kernel parameter
#---------------------------------------------------
#' Default mutator of the field \code{Kernpar} of the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param kernpar [\code{vector(1)}]\cr
#' New kernel parameters.
#'
#' @return The object taken in the argument.
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setKernpar.default <- function(kern, kernpar) {
  return(kern)
}

#--- Specific mutator method for the kernel parameter
#----------------------------------------------------
#' Mutator of the field \code{Kernpar} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel} taken in the argument.
#' @param kernpar [\code{vector(1)}]\cr
#' New kernel parameters.
#'
#' @return [\code{\link{Kernel}(1)}]
#' Object of class \code{Kernel} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab Type of kernel, \cr
#'    \code{Mat} \tab Kernel matrix, \cr
#'    \code{Kernpar} \tab Parameters of kernel, when required, \cr
#'    \code{bincat} \tab Index of binary/categorical variables, when required,
#'  }
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setKernpar.Kernel<-function(kern, kernpar) {
  kern$Kernpar <- kernpar
  return(kern)
}

#--- Generic mutator method for the kernel parameter bincat: just concerning the additive kernel
#---------------------------------------------------
#' Default mutator of the field \code{bincat} of the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel} taken in the argument.
#' @param bincat [\code{vector(1)}]\cr
#' Index of binary/categorial variables.
#'
#' @return The object taken in the argument.
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setBincat <- function(kern, bincat) {
  UseMethod("setMat", bincat)
}
#--- Default mutator method for the kernel parameter bincat: just concerning the additive kernel
#-----------------------------------------------------------------------------------------------
#' Default mutator of the field \code{bincat} of the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @param bincat [\code{vector(1)}]\cr
#' New index of binary/categorial variables.
#'
#' @return The object taken in the argument.
#'
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
setBincat.default <- function(kern, bincat) {
  return(bincat)
}

#--- Specific mutator method for the kernel parameter bincat: just concerning the additive kernel
#----------------------------------------------------
#' Mutator of the field \code{bincat} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel} taken in the argument.
#' @param bincat [\code{vector(1)}]\cr
#' New index of binary/categorial variables.
#'
#' @return [\code{\link{Kernel}(1)}]
#' Object of class \code{Kernel} with elements:
#'  \tabular{ll}{
#'    \code{Type} \tab type of kernel, \cr
#'    \code{Mat} \tab kernel matrix, \cr
#'    \code{Kernpar} \tab parameters of kernel, when required, \cr
#'    \code{bincat} \tab index of binary/categorical variables, when required.
#'  }
#'
#' @author Cesaire J. K. Fouodo
#'
#' @keywords internal
setBincat.Kernel<-function(kern, bincat) {
  kern$bincat <- bincat
  return(kern)
}


#>>>> Accessors for the kernel object

#--- Generic accessor method for the kernel type
#-----------------------------------------------
#' Creator of the generic accessor \code{getType}.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getType <- function(kern) {
  UseMethod("getType", kern)
}

#--- Default accessor method for the kernel type
#-----------------------------------------------
#' Accessor for the field \code{Type} for the object taken in the argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getType.default <- function(kern) {
  return(NULL)
}

#--- Specific accessor method for the kernel type
#------------------------------------------------
#' Accessor for the field \code{Type} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#'
#' @return Type of the kernel taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getType.Kernel <- function(kern) {
  return(kern$Type)
}

#--- Generic accessor method for the kernel matrix
#-----------------------------------------------
#' Creator of the generic accessor \code{getMat}.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getMat <- function(kern) {
  UseMethod("getMat", kern)
}
#--- Default accessor method for the kernel matrix
#-----------------------------------------------
#' Accessor for the field \code{Mat} for the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#'
#' @author Cesaire J. K. Fouodo
getMat.default <- function(kern) {
  return(NULL)
}

#--- Specific accessor method for the kernel matrix
#------------------------------------------------
#' Accessor for the field \code{Mat} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#'
#' @return [\code{matrix(1)}]
#' The kernel matrix.
#'
#' @keywords internal
#' @author Cesaire J. K. Fouodo
getMat.Kernel<-function(kern) {
  return(kern$Mat)
}

#--- Generic accessor method for the kernel parameter
#---------------------------------------------------
#' Creator of the generic accessor \code{getKernpar}.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernpar <- function(kern) {
  UseMethod("getKernpar", kern)
}
#--- Default accessor method for the kernel parameter
#-----------------------------------------------
#' Accessor for the field \code{Kernpar} for the object taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernpar.default <- function(kern) {
  return(NULL)
}

#--- Specific accessor method for the kernel parameter
#------------------------------------------------
#' Accessor for the field \code{Kernpar} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#'
#' @return [\code{vector(1)}]
#' The kernel parameters.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getKernpar.Kernel <- function(kern) {
  return(kern$Kernpar)
}

#--- Generic accessor method for the kernel bincat parameter: just concerning the additive kernel
#-----------------------------------------------
#' Creator of the generic accessor \code{getBincat}.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBincat <- function(kern) {
  UseMethod("getKernpar", kern)
}
#--- Default accessor method for the kernel parameter: just concerning the additive kernel
#-----------------------------------------------
#' Accessor for the field \code{Bincat} for the object taken in the argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object taken in the argument.
#'
#' @return \code{NULL}.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBinca.default <- function(kern) {
  return(NULL)
}

#--- Specific accessor method for the kernel parameter: just concerning the additive kernel
#------------------------------------------------
#' Accessor for the field \code{Bincat} of the object of class \code{Kernel} taken in an argument.
#'
#'
#' @title \code{Kernel}
#' @param kern [\code{\link{Kernel}(1)}]\cr
#' Object of class \code{Kernel}.
#'
#' @return Index of binary/categorical variables.
#' @keywords internal
#'
#' @author Cesaire J. K. Fouodo
getBincat.Kernel <- function(kern) {
  return(kern$bincat)
}
