\name{orlm}
\alias{orlm}
\alias{orlm.formula}

\title{Fitting multivariate regression models with order restrictions}
\description{It determines the order-restricted maximum likelihood estimates 
and the corresponding log likelihood for the hypothesis of interest.
Additionally it gives the (unconstrained) maximum likelihood estimates and
the active contraints.}

\usage{
orlm(formula, data, constr, rhs, nec, control = orlmcontrol())
\method{orlm}{formula}(formula, data, constr, rhs, nec, control = orlmcontrol())
}
\arguments{
  \item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
}
  \item{constr}{matrix with constraints; with rows as constraint
    definition, columns should be in line with the parameters of the model}
  \item{rhs}{vector of right hand side elements; \eqn{Constr \; \theta \geq
      rhs}; number should equal the number of rows of the constr matrix}
  \item{nec}{number of equality constraints; all remaining ones are
    inequality constraints; should be \eqn{\leq} the length of rhs}
  \item{control}{a list of control arguments; see \code{\link{orlmcontrol}} for details.
  }
  }
\details{
The contraints in the hypothesis of interest are defined by 
\eqn{Constr}, \eqn{rhs}, and \eqn{nec}.
The first \eqn{nec} constraints are the equality contraints:
\eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]};
and the remaing ones are the inequality contraints:
\eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}.

Two requirements should be met:
\enumerate{
\item The first \eqn{nec} constraints must be the equality contraints 
(i.e., \eqn{Constr[1:nec, 1:tk] \theta = rhs[1:nec]}) 
and the remaining ones the inequality contraints
(i.e., \eqn{Constr[nec+1:c_m, 1:tk] \theta \geq rhs[nec+1:c_m]}).
\item When \eqn{rhs} is not zero, \eqn{Constr} should be of full rank 
(after discarding redundant restrictions). 
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011) and Kuiper, Hoijtink, and Silvapulle (unpublished).
}
}
\value{an object of class orlm}
\references{
Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (unpublished). Generalization of the 
Order-Restricted Information Criterion for Multivariate Normal Linear Models.
}
\author{Daniel Gerhard and Rebecca M. Kuiper}


\seealso{\code{\link{solve.QP}}, \code{\link{goric}}}

\examples{
## Example from Kuiper, R.M. and Hoijtink, H. (Unpublished). A Fortran 90 program for the generalization of the order restricted information criterion.

# constraint definition
cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
constr <- kronecker(diag(3), cmat)

# no effect model
fm0 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene, constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr))
fm0

# order constrained model (increasing serum levels with increasing doses)
fm1 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene, constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr))
summary(fm1)

# unconstrained model
fmunc <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene, constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0)
fmunc
}

\keyword{models}
\keyword{optimize}
