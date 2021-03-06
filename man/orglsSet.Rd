% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orglsSet.R
\name{orglsSet}
\alias{orglsSet}
\title{Set of generalised least-squares models}
\usage{
orglsSet(formula, data, weights = NULL, correlation = NULL, set,
  direction = "increase", n = NULL, base = 1, control = orlmcontrol())
}
\arguments{
\item{formula}{an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.}

\item{weights}{a \code{\link{varClasses}} object; more details are provided on the help pages in R package \code{nlme}}

\item{correlation}{a \code{\link{corClasses}} object; more details are provided on the help pages in R package \code{nlme}}

\item{set}{either a character string (see \code{\link{constrSet}}), or a list with slots for constr, rhs, and nec similarly defined as in \code{\link{orlm}}}

\item{direction}{direction of the order constraints}

\item{n}{a (possibly named) vector of sample sizes for each group}

\item{base}{column of the constraint matrix representing a control group}

\item{control}{a list of control arguments; see \code{\link{orlmcontrol}} for details.}
}
\value{
a list with orgls objects
}
\description{
Fitting a specific set of generalisd least-squares models with order restrictions.
}
\details{
This function is just a wrapper for repeated calls of \code{\link{orgls}} with different constraint definitions. Predefined lists with constraint-sets can be constructed with function \code{\link{constrSet}}.
}
\seealso{
\code{\link{orgls}}, \code{\link{constrSet}}, \code{\link{goric}}
}
\keyword{models}
