#' Performs support vectors analysis for data sets with survival outcome. 
#' Three approaches are available in the package: 
#' The regression approach takes censoring into account when formulating the inequality constraints of the support vector problem. 
#' In the ranking approach, the inequality constraints set the objective to maximize the concordance index for comparable pairs of observations. 
#' The hybrid approach combines the regression and ranking constraints in the same model.
#' 
#' The following denotations are used for the models implemented:
#'  \itemize{
#'    \item \code{'regression'} referring to the regression approach, named \code{SVCR} in Van Belle et al. (2011b),
#'    \item \code{'vanbelle1'} according to the first version of survival surpport vector machines based on ranking constraints,
#'          named \code{RANKSVMC} by Van Belle et al. (2011b),
#'    \item \code{'vanbelle2'} according to the second version of survival surpport vector machines based on ranking constraints
#'          like presented in \code{model1} by Van Belle et al. (2011b) and
#'    \item \code{'hybrid'}  combines simultaneously the regression and ranking constraints in the same model. Hybrid model is labeled
#'          \code{model2} by Van Belle et al. (2011b).
#'
#'  }
#' The argument \code{'type'} of the function \code{survivalsvm} is used to set the type of model to be fitted.
#' For the models \code{vanbelle1}, \code{vanbelle2} and \code{hybrid}, differences between comparable
#' pairs of observations are required. Each observation is compared with its nearest neighbor according to the survival time, and the
#' three possible comparison approaches \link{makediff1}, \link{makediff2} and \link{makediff3} are offered to compute the
#' differences between comparable neighbors.
#'
#' The current version of \code{survivalsvm} uses the solvers \code{\link{ipop}} and \code{\link{quadprog}} to solve the dual
#' optimization problems deduced from the suport vector formulations of the models presented above. Notice that for using \code{quadprog}
#' the kernel matrix needs to be symmetric and positive definite. Therefore when the conditions are not met, the kernel matrix needs be slightly perturbed to obtain the nearest positive definite kernel matrix.
#' The alternative to \code{quadprog} is \code{ipop}, that can also handle a non-negative definite kernel matrix, however more time may be
#' required to solve the quadratic optimization dual problem. The argument \code{opt.meth} is used to select the solver.
#'
#' The \code{survivalsvm} command can be called giving a formula, in which the survival time and the status are grouped into a
#' two colunm matrix using the command \code{\link{Surv}} from the package \code{survival}. An alternative is to pass the data
#' frame of training data points as an argument using \code{data}, to mention the name of the survival time variable and
#' the name of the status variable as illustrated in the third example below.
#'
#' @title survivalsvm
#' @param formula [\code{formula(1)}]\cr
#' Object of class \code{formula}. See \code{\link{formula}} for more details.
#' @param data [\code{data.frame(1)}]\cr
#' Object of class \code{data.frame} containing data points that will be used to fit the model.
#' @param subset [\code{vector(1)}]\cr
#' An index vector specifying the cases to be used in the training sample.
#' @param type [\code{character(1)}]\cr
#' String indicating which type of survival support vectors model is desired. This must be one
#' of the following strings: 'regression', 'vanbelle1', 'vanbelle2' or 'hybrid'.
#' @param diff.meth [\code{character(1)}]\cr
#' String indicating which of \code{'makediff1'}, \code{'makediff2'} or \code{'makediff3'}
#' is used in case of 'vanbelle1', 'vanbelle2' and 'hybrid'.
#' @param gamma.mu [\code{numeric(1)|vector(1)}]\cr
#' Parameters of regularization. Note that a vector with two parameters is required in case of \code{hybrid} approach. Just
#' one value is required in case of \code{regression}, \code{vanbelle1} or \code{vanbelle2}.
#' @param opt.meth [\code{character(1)}]\cr
#' Program used to solve the quadratic optimization problem. Either "\code{\link{quadprog}}" or "\code{\link{ipop}}".
#' @param kernel [\code{\link{Kernel}(1)}]\cr
#' Kernel used to fit the model: linear kern ('lin_kernel'), additive kernel ('add_kernel'),
#' radial basis kernels ('rbf_kernel') and the polynomial kernel ('poly_kernel').
#' @param kernel.pars [\code{vector(1)}]\cr
#' Parameters of kernel, when required.
#' @param time.variable.name [\code{character}]\cr
#' Name of the survival time variable in \code{data}, when given in argument.
#' @param status.variable.name [\code{character(1)}]\cr
#' Name of the status variable in \code{data}.
#' @param sgf.sv [\code{character(1)}]\cr
#' Number of decimal digits in the solution of the quadratic optimization problem.
#' @param sigf [\code{numeric(1)}]\cr
#' Used by \code{\link{ipop}}. See \code{\link{ipop}} for details.
#' @param maxiter [\code{integer(1)}]\cr
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
#' @return \code{survivalsvm}
#' Object of class \code{survivalsvm}, with elements:
#' \tabular{ll}{
#'    \code{call} \tab command calling this program, \cr
#'    \code{typeofsurvivalsvm} \tab type of survival support vector machines approach,\cr
#'    \code{model.fit} \tab the fitted survival model,\cr
#'    \code{var.names} \tab names of variables used.\cr
#'  }
#' @export
#'
#' @seealso \link{predict.survivalsvm}
#' @examples 
#' survivalsvm(Surv(time, status) ~ ., veteran, gamma.mu = 0.1)
#'
#' survsvm.reg <- survivalsvm(formula = Surv(diagtime, status) ~ ., data = veteran,
#'                            type = "regression", gamma.mu = 0.1,
#'                            opt.meth = "ipop", kernel = "add_kernel")
#'                             
#' survsvm.vb2 <- survivalsvm(data = veteran, time.variable.name = "diagtime",
#'                            status.variable.name = "status", 
#'                            type = "vanbelle2", gamma.mu = 0.1,
#'                            opt.meth = "quadprog", diff.meth = "makediff3", 
#'                            kernel = "lin_kernel",
#'                            sgf.sv = 5, sigf = 7, maxiter = 20, 
#'                            margin = 0.05, bound = 10)
#'                             
#' @author Cesaire J. K. Fouodo
#'
#' @note This implementation is in part inspired by the \code{Matlab} toolbox \code{Survlab}
#'  (\href{http://user.it.uu.se/~kripe367/survlab/instruction.html}{\code{A Survival Analysis Toolbox}}).
#' @references
#' \itemize{
#'  \item Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011a).
#'        Improved performance on high-dimensional survival data by application of Survival-SVM.
#'        Bioinformatics (Oxford, England) 27, 87-94.
#'  \item Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011b).
#'        Support vector methods for survival analysis: a comparaison between ranking and regression approaches.
#'        Artificial Intelligence in medecine 53, 107-118.
#' }
#' @import survival
#' @importFrom utils packageVersion
survivalsvm <- function (formula = NULL, data = NULL, subset = NULL, type = "regression",
                        diff.meth = NULL, gamma.mu = NULL, opt.meth = "quadprog", kernel = "lin_kernel",
                        kernel.pars = NULL, time.variable.name = NULL, status.variable.name = NULL,
                        sgf.sv = 5, sigf = 7, maxiter = 20, margin = 0.05, bound = 10,
                        eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08) {
  if (!(tolower(type) %in% c("regression", "vanbelle1", "vanbelle2", "hybrid"))) {
    stop("Error: 'type' must be either 'regression', 'vanbelle1', 'vanbelle2' or 'hybrid'.")
  }
  if (!is.null(diff.meth)) {
    if (tolower(type) != "regression") {
      if (!(tolower(diff.meth) %in% c("makediff1", "makediff2", "makediff3"))) {
        stop("'diff.meth' must be either 'makediff1', 'makediff2' or 'makediff3'.")
      }
    }
  } else {
    if (tolower(type) %in% c("vanbelle1", "vanbelle2", "hybrid")) {
      stop("Types 'vanbelle1', 'vanbell2' and 'hybrid' require an argument diff.meth.")
    }
  }
  if (is.null(gamma.mu)) {
    stop("gamma.mu can not be NULL.")
  }
  if (any(gamma.mu <= 0)) {
    stop("gamma.mu: only positive values allowed.")
  }
  if (tolower(type) == "hybrid") {
    if(length(gamma.mu) != 2) {
      stop("'gamma.mu' must be a vector of two numeric values.")
    }
  } else {
    if(length(gamma.mu) == 2) {
      warning("The second element of 'gamma.mu' has been ignored because of not hybrid type.")
    }
  }
  if (is.null(kernel.pars)){
    kernel.pars <- NA
  }
  if (!(tolower(opt.meth) %in% c("quadprog", "ipop"))) {
    stop("Error: opt.meth must be  either 'quadprog' or 'ipop'.")
  }
  if (!(tolower(kernel) %in% c("lin_kernel", "add_kernel", "rbf_kernel", "rbf4_kernel", "poly_kernel"))) {
    stop("'kernel' must be either 'lin_kern', 'add_kernel', 'rbf_kernel', 'rbf4_kernel' or 'poly_kernel'.")
  }
  if (!inherits(sgf.sv, "numeric")) {
    stop("'sgf.sv' must be a numeric.")
  } else {
    if(sgf.sv <= 0) {
      stop("'sgf.sv' must be greather than 0.")
    }
  }
  # test on the formula
  if (is.null(formula)) {
    if (is.null(time.variable.name)) {
      stop("Error: Please give a formula or dependent variable names.")
    }
    if (is.null(data)) {
      stop("'data' can not be NULL.")
    }
    if (!(inherits(data, "data.frame"))) {
      stop("data must be a data.frame")
    }
    if (is.null(status.variable.name)) {
      status.variable.name <- "none"
      response <- data[, time.variable.name]
    } else {
      if (is.null(subset)) {
        # the given data is a data frame
        response <- data[, c(time.variable.name, status.variable.name)]
        covar <- setdiff(names(data), c(time.variable.name, status.variable.name))
        if (length(covar) == 0)
          stop("No covariable found.")
        traindata <- data[, setdiff(names(data), c(time.variable.name, status.variable.name))]
      } else {
        response <- data[subset, c(time.variable.name, status.variable.name)]
        covar <- setdiff(names(data), c(time.variable.name, status.variable.name))
        if (length(covar) == 0) {
          stop("No covariable hat been found.")
        }
        traindata <- data[subset, setdiff(names(data), c(time.variable.name, status.variable.name))]
      }
      index.factor <- if (kernel == "add_kernel") which(sapply(traindata, is.factor)) else integer(0)
      X <- data.matrix(traindata)
      X <- X[stats::complete.cases(X), , drop = FALSE]
      if (prod(dim(X)) == 0) {
        stop("No observation in the data frame given in argument.")
      }
      Y <- if (status.variable.name == "none") response else response[,1]
      delta <- if (status.variable.name == "none") rep(1, length(Y)) else response[, 2]
    }
  } else {
    formula <- formula(formula)
    if (!inherits(formula,"formula")) {
      stop("Error: Invalid formula.")
    }
    data.selected <- stats::model.frame(formula, data)
    response <- data.selected[[1]]
    if (inherits(response, "matrix")) {
      is.named <- FALSE
      if (!is.null(status.variable.name) && !is.null(time.variable.name)) {
        if (status.variable.name %in% names(response)) {
          delta <- response[, status.variable.name]
          if (time.variable.name %in% names(response)) {
            Y <- response[, time.variable.name]
            is.named <- TRUE
          }
        }
      }
      if (!is.named)  {
      if (ncol(response)  > 2) {
        stop("Error: Please names the intresting column by 'time' and 'status'.")
        }
      Y <- response[, 1]
      delta <- response[, 2]
      }
    }# the response was a matrix
    if (inherits(response, "Surv")) {
      response <- as.matrix(response)
        Y <- response[, "time"]
        delta <- response[, "status"]
      }
    traindata <- data.selected[-1]
    # extract the subset of interest data points, when given
    if (!is.null(subset)) {
      if (!all(subset %in% 1:length(Y))) {
        stop("Error: invalid subset.")
      }
      Y <- Y[subset]
      delta <- delta[subset]
      traindata <- traindata[subset, , drop = FALSE]
    }
    covar <- names(traindata)
    index.factor <- if(kernel == "add_kernel") which(sapply(traindata, is.factor)) else integer(0)
    X <- data.matrix(traindata)
    X <- X[stats::complete.cases(X), , drop = FALSE]
    if (ncol(X) == 0) {
      stop("No variable in the data frame given in argument.")
    }
    if (nrow(X) == 0) {
      stop("No observation in the data frame given in argument.")
    }
  }
  if (!inherits(Y, "numeric")) {
    stop("The time must be a numeric vector.")
  }
  if (!all(delta %in% c(0,1))) {
    stop("Error: Status must either be 0 or 1.")
  }
  #selection of a method
  if (tolower(type) == "regression") {
    model.fit <- regFit(X = X, Y = Y, delta = delta, meth_par = gamma.mu, kernel_type = kernel,
                        kernel_pars = kernel.pars, bin_cat = index.factor, opt_alg = opt.meth,
                        sgf_sv = sgf.sv, sigf = sigf, maxiter = maxiter, margin = margin, bound = bound,
                        eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)
  }
  if (tolower(type) == "vanbelle1") {
    model.fit <- vanbelle1Fit(X = X, Y = Y, delta = delta,
                              meth_par = gamma.mu, kernel_type = kernel,
                              kernel_pars = kernel.pars, bin_cat = index.factor,
                              makediff = match.fun(diff.meth), opt_alg = opt.meth, sgf_sv = sgf.sv,
                              sigf = sigf, maxiter = maxiter, margin = margin, bound = bound,
                              eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)
  }
  if (tolower(type) == "vanbelle2") {
    model.fit <- vanbelle2Fit(X = X, Y = Y, delta = delta,
                              meth_par = gamma.mu, kernel_type = kernel,
                              kernel_pars = kernel.pars, bin_cat = index.factor,
                              makediff = match.fun(diff.meth), opt_alg = opt.meth, sgf_sv = sgf.sv,
                              sigf = sigf, maxiter = maxiter, margin = margin, bound = bound,
                              eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)
  }
  if (tolower(type) == "hybrid") {
    model.fit <- hybridFit(X = X, Y = Y, delta = delta,
                           meth_par = gamma.mu, kernel_type = kernel,
                           kernel_pars = kernel.pars, bin_cat = index.factor,
                           makediff = match.fun(diff.meth), opt_alg = opt.meth, sgf_sv = sgf.sv,
                           sigf = sigf, maxiter = maxiter, margin = margin, bound = bound,
                           eig.tol = eig.tol, conv.tol = conv.tol, posd.tol = posd.tol)
  }
  result <- list(
    call = sys.call(),
    typeofsurvivalsvm = tolower(type),
    model.fit = model.fit,
    var.names = covar
  )
  class(result) <- "survivalsvm"
  result$package.version <- unlist(packageVersion("survivalsvm"))
  return(result)
}
