% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survivalsvm.R
\name{survivalsvm}
\alias{survivalsvm}
\title{survivalsvm}
\usage{
survivalsvm(
  formula = NULL,
  data = NULL,
  subset = NULL,
  type = "regression",
  diff.meth = NULL,
  gamma.mu = NULL,
  opt.meth = "quadprog",
  kernel = "lin_kernel",
  kernel.pars = NULL,
  time.variable.name = NULL,
  status.variable.name = NULL,
  sgf.sv = 5,
  sigf = 7,
  maxiter = 20,
  margin = 0.05,
  bound = 10,
  eig.tol = 1e-06,
  conv.tol = 1e-07,
  posd.tol = 1e-08
)
}
\arguments{
\item{formula}{[\code{formula(1)}]\cr
Object of class \code{formula}. See \code{\link{formula}} for more details.}

\item{data}{[\code{data.frame(1)}]\cr
Object of class \code{data.frame} containing data points that will be used to fit the model.}

\item{subset}{[\code{vector(1)}]\cr
An index vector specifying the cases to be used in the training sample.}

\item{type}{[\code{character(1)}]\cr
String indicating which type of survival support vectors model is desired. This must be one
of the following strings: 'regression', 'vanbelle1', 'vanbelle2' or 'hybrid'.}

\item{diff.meth}{[\code{character(1)}]\cr
String indicating which of \code{'makediff1'}, \code{'makediff2'} or \code{'makediff3'}
is used in case of 'vanbelle1', 'vanbelle2' and 'hybrid'.}

\item{gamma.mu}{[\code{numeric(1)|vector(1)}]\cr
Parameters of regularization. Note that a vector with two parameters is required in case of \code{hybrid} approach. Just
one value is required in case of \code{regression}, \code{vanbelle1} or \code{vanbelle2}.}

\item{opt.meth}{[\code{character(1)}]\cr
Program used to solve the quadratic optimization problem. Either "\code{\link[pracma]{quadprog}}" or "\code{\link[kernlab]{ipop}}".}

\item{kernel}{[\code{\link{Kernel}(1)}]\cr
Kernel used to fit the model: linear kern ('lin_kernel'), additive kernel ('add_kernel'),
radial basis kernels ('rbf_kernel') and the polynomial kernel ('poly_kernel').}

\item{kernel.pars}{[\code{vector(1)}]\cr
Parameters of kernel, when required.}

\item{time.variable.name}{[\code{character}]\cr
Name of the survival time variable in \code{data}, when given in argument.}

\item{status.variable.name}{[\code{character(1)}]\cr
Name of the status variable in \code{data}.}

\item{sgf.sv}{[\code{character(1)}]\cr
Number of decimal digits in the solution of the quadratic optimization problem.}

\item{sigf}{[\code{numeric(1)}]\cr
Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.}

\item{maxiter}{[\code{integer(1)}]\cr
Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.}

\item{margin}{[\code{numeric(1)}]\cr
Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.}

\item{bound}{[\code{numeric(1)}]\cr
Used by \code{\link[kernlab]{ipop}}. See \code{\link[kernlab]{ipop}} for details.}

\item{eig.tol}{[\code{numeric(1)}]\cr
Used by \code{\link[Matrix]{nearPD}} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.}

\item{conv.tol}{[\code{numeric(1)}]\cr
Used by \code{\link[Matrix]{nearPD}} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.}

\item{posd.tol}{[\code{numeric(1)}]\cr
Used by \code{\link[Matrix]{nearPD}} for adjusting positive definiteness. See \code{\link[Matrix]{nearPD}} for detail.}
}
\value{
\code{survivalsvm}
Object of class \code{survivalsvm}, with elements:
\tabular{ll}{
   \code{call} \tab command calling this program, \cr
   \code{typeofsurvivalsvm} \tab type of survival support vector machines approach,\cr
   \code{model.fit} \tab the fitted survival model,\cr
   \code{var.names} \tab names of variables used.\cr
 }
}
\description{
Performs support vectors analysis for data sets with survival outcome. 
Three approaches are available in the package: 
The regression approach takes censoring into account when formulating the inequality constraints of the support vector problem. 
In the ranking approach, the inequality constraints set the objective to maximize the concordance index for comparable pairs of observations. 
The hybrid approach combines the regression and ranking constraints in the same model.
}
\details{
The following denotations are used for the models implemented:
 \itemize{
   \item \code{'regression'} referring to the regression approach, named \code{SVCR} in Van Belle et al. (2011b),
   \item \code{'vanbelle1'} according to the first version of survival surpport vector machines based on ranking constraints,
         named \code{RANKSVMC} by Van Belle et al. (2011b),
   \item \code{'vanbelle2'} according to the second version of survival surpport vector machines based on ranking constraints
         like presented in \code{model1} by Van Belle et al. (2011b) and
   \item \code{'hybrid'}  combines simultaneously the regression and ranking constraints in the same model. Hybrid model is labeled
         \code{model2} by Van Belle et al. (2011b).

 }
The argument \code{'type'} of the function \code{survivalsvm} is used to set the type of model to be fitted.
For the models \code{vanbelle1}, \code{vanbelle2} and \code{hybrid}, differences between comparable
pairs of observations are required. Each observation is compared with its nearest neighbor according to the survival time, and the
three possible comparison approaches \link{makediff1}, \link{makediff2} and \link{makediff3} are offered to compute the
differences between comparable neighbors.

The current version of \code{survivalsvm} uses the solvers \code{\link[kernlab]{ipop}} and \code{\link[pracma]{quadprog}} to solve the dual
optimization problems deduced from the suport vector formulations of the models presented above. Notice that for using \code{quadprog}
the kernel matrix needs to be symmetric and positive definite. Therefore when the conditions are not met, the kernel matrix needs be slightly perturbed to obtain the nearest positive definite kernel matrix.
The alternative to \code{quadprog} is \code{ipop}, that can also handle a non-negative definite kernel matrix, however more time may be
required to solve the quadratic optimization dual problem. The argument \code{opt.meth} is used to select the solver.

The \code{survivalsvm} command can be called giving a formula, in which the survival time and the status are grouped into a
two colunm matrix using the command \code{\link[survival]{Surv}} from the package \code{survival}. An alternative is to pass the data
frame of training data points as an argument using \code{data}, to mention the name of the survival time variable and
the name of the status variable as illustrated in the third example below.
}
\note{
This implementation is in part inspired by the \code{Matlab} toolbox \code{Survlab}
 (\href{http://user.it.uu.se/~kripe367/survlab/instruction.html}{\code{A Survival Analysis Toolbox}}).
}
\examples{
survivalsvm(Surv(time, status) ~ ., veteran, gamma.mu = 0.1)

survsvm.reg <- survivalsvm(formula = Surv(diagtime, status) ~ ., data = veteran,
                           type = "regression", gamma.mu = 0.1,
                           opt.meth = "ipop", kernel = "add_kernel")
                            
survsvm.vb2 <- survivalsvm(data = veteran, time.variable.name = "diagtime",
                           status.variable.name = "status", 
                           type = "vanbelle2", gamma.mu = 0.1,
                           opt.meth = "quadprog", diff.meth = "makediff3", 
                           kernel = "lin_kernel",
                           sgf.sv = 5, sigf = 7, maxiter = 20, 
                           margin = 0.05, bound = 10)
                            
}
\references{
\itemize{
 \item Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011a).
       Improved performance on high-dimensional survival data by application of Survival-SVM.
       Bioinformatics (Oxford, England) 27, 87-94.
 \item Van Belle, V., Pelcmans, K., Van Huffel S. and Suykens J. A.K. (2011b).
       Support vector methods for survival analysis: a comparaison between ranking and regression approaches.
       Artificial Intelligence in medecine 53, 107-118.
}
}
\seealso{
\link{predict.survivalsvm}
}
\author{
Cesaire J. K. Fouodo
}
