% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\docType{class}
\name{model-class}
\alias{model-class}
\alias{.model}
\title{Finmix \code{model} class}
\description{
This class specifies a finite mixture model. Entities are created from it by
calling its constructor \code{\link[=model]{model()}}.
}
\details{
A finite mixture model in the \code{ finmix} package is defined by its number of
components \code{K}, the component distributions \code{dist}, the data dimension \code{r}
and an indicator defining, if the model has fixed or unknown indicators.
Finite mixture models for the following distributions can be constructed:
\itemize{
\item Poisson,
\item Conditional Poisson,
\item Exponential,
\item Binomial,
\item Normal,
\item Multivariate Normal,
\item Student-t,
\item Multivariate Student-t.
}

Using the constructor \code{\link[=model]{model()}} a finite mixture model can be created, the
default being a mixture model of Poisson distributions.
\subsection{Fully defined finite mixture models}{

A fully defined finite mixture model contains next to the distribution and
the components also weights and parameters. The weights are defined in slot
\code{weight} and must be of class \code{ matrix} with as many weights as there are
components in the mixture model (dimension \code{Kx1}). Parameters are defined in
a \code{ list} named \code{par}. The elements of this list depend on the chosen
distribution in slot \code{dist}:
\itemize{
\item Poisson: A \code{matrix} named \code{lambda} of dimension \code{Kx1} holding the rate
parameters.
\item Exponential: A \code{matrix} named \code{lambda} of dimension \code{Kx1} holding the rate
parameters.
\item Binomial: A \code{matrix} of dimension \code{Kx1} named \code{p} storing the
probabilities.
}
}
}
\section{Slots}{

\describe{
\item{\code{dist}}{A character, defining the distribution family. Possible choices
are binomial, exponential, normal, normult, poisson, student, and studmult.}

\item{\code{r}}{An integer. Defines the vector dimension of a model. Is one for all
univariate distributions and larger than one for normult and studmult.}

\item{\code{K}}{An integer, defining the number of components in the finite mixture.}

\item{\code{weight}}{A matrix, containing the weights of the finite mixture model.
The matrix must have dimension \code{1 x K} and weights must add to one
must all be larger or equal to zero.}

\item{\code{par}}{A list containing the parameter vectors for the finite mixture
distribution. The list can contain more than one named parameter vector.}

\item{\code{indicmod}}{A character defining the indicator model. So far only
multinomial indicator models are possible.}

\item{\code{indicfix}}{A logical. If \code{TRUE} the indicators are given and
therefore fixed.}

\item{\code{T}}{A matrix containing the repetitions in case of a \code{"binomial"} or
\code{"poisson"} model.}
}}

\seealso{
\itemize{
\item \code{\link[=mixturemcmc]{mixturemcmc()}} for performing MCMC sampling with a mixture model
\item \code{\link[=modelmoments]{modelmoments()}} for compute theoretical moments of a finite mixture model
}
}
