\name{goric}
\alias{goric}
\alias{goric_penalty}

\title{Calculate GORIC}
\description{
It calculates the order-restricted log likelihood, the penalty of the GORIC,
the GORIC values, and the GORIC weights for a set of hypotheses (orlm objects),
where the penalty is based on \eqn{iter} iterations.
}
\usage{
goric(object, ..., iter = 1e+05, mc.cores=1)
goric_penalty(object, iter = 1e+05, mc.cores=1)
}
\arguments{
  \item{object}{an object of class orlm}
  \item{\dots}{further objects of class orlm}
  \item{iter}{number of iterations to calculate GORIC penalty terms}
  \item{mc.cores}{number of cores using a socket cluster implemented in package
  \code{parallel}}
}
\details{
The hypothesis with the lowest GORIC value is the preferred one.

The GORIC weights reflect the support of each hypothesis in the set.
To compare two hypotheses (and not one to the whole set), one should examine 
the ratio of the two corresponding GORIC weights.


To safequard for weak hypotheses (i.e., hypotheses not supported by the data), 
one should include a model with no constraints (the so-called unconstrained model).
More information can be obtained from Kuiper, Hoijtink, and Silvapulle (2011) and 
Kuiper, Hoijtink, and Silvapulle (unpublished).
}
\value{a data.frame with the information criteria or a single penalty term}
\references{
Kuiper R.M., Hoijtink H., Silvapulle M.J. (2011). An Akaike-type Information
Criterion for Model Selection Under Inequality Constraints. \emph{Biometrika}, 
\bold{98}, 495--501.

Kuiper R.M., Hoijtink H., Silvapulle M.J. (unpublished). Generalization of the 
Order-Restricted Information Criterion for Multivariate Normal Linear Models.
}
\author{
Daniel Gerhard and Rebecca M. Kuiper
}


\seealso{\code{\link{orlm}}}
\examples{
## Example from Kuiper, R.M. and Hoijtink, H. (Unpublished).
# A Fortran 90 program for the generalization of the
# order restricted information criterion.

# constraint definition
cmat <- cbind(diag(3), 0) + cbind(0, -diag(3))
constr <- kronecker(diag(3), cmat)
constr

# no effect model
fm0 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=nrow(constr))
fm0

# order constrained model (increasing serum levels with increasing doses)
fm1 <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
            constr=constr, rhs=rep(0, nrow(constr)), nec=0)
summary(fm1)

# unconstrained model
fmunc <- orlm(cbind(SDH, SGOT, SGPT) ~ dose-1, data=vinylidene,
              constr=matrix(0, nrow=1, ncol=12), rhs=0, nec=0)
fmunc

# calculate GORIC
# (only small number of iterations to decrease computation time, default: iter=100000)
goric(fm0, fm1, fmunc, iter=1000)
}

\keyword{models}
\keyword{optimize}
