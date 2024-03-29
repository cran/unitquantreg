% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpqr-uburrxii.R
\name{uburrxii}
\alias{uburrxii}
\alias{duburrxii}
\alias{puburrxii}
\alias{quburrxii}
\alias{ruburrxii}
\title{The unit-Burr-XII distribution}
\usage{
duburrxii(x, mu, theta, tau = 0.5, log = FALSE)

puburrxii(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

quburrxii(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

ruburrxii(n, mu, theta, tau = 0.5)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{mu}{location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.}

\item{theta}{nonnegative shape parameter.}

\item{tau}{the parameter to specify which quantile is to used.}

\item{log, log.p}{logical; If TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\value{
\code{duburrxii} gives the density, \code{puburrxii} gives the distribution function,
\code{quburrxii} gives the quantile function and \code{ruburrxii} generates random deviates.

Invalid arguments will return an error message.
}
\description{
Density function, distribution function, quantile function and random number generation function
for the unit-Burr-XII distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
}
\details{
Probability density function
\deqn{f(y\mid \alpha, \theta )=\frac{\alpha \theta }{y}\left[ -\log (y)\right]^{\theta -1}\left\{ 1+\left[ -\log (y)\right] ^{\theta }\right\} ^{-\alpha -1}}

Cumulative distribution function
\deqn{F(y\mid \alpha, \theta )=\left\{ 1+\left[ -\log (y)\right] ^{\theta}\right\} ^{-\alpha }}

Quantile function
\deqn{Q(\tau \mid \alpha, \theta )=\exp \left[ -\left( \tau ^{-\frac{1}{\alpha }}-1\right)^{\frac{1}{\theta }} \right]}

Reparameterization
\deqn{\alpha=g^{-1}(\mu)=\frac{\log\left ( \tau^{-1} \right )}{\log\left [ 1+\log\left ( \frac{1}{\mu} \right )^\theta \right ]}}
}
\examples{
set.seed(123)
x <- ruburrxii(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
R <- range(x)
S <- seq(from = R[1], to = R[2], by =  0.01)
hist(x, prob = TRUE, main = 'unit-Burr-XII')
lines(S, duburrxii(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(ecdf(x))
lines(S, puburrxii(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(quantile(x, probs = S), type = "l")
lines(quburrxii(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)

}
\references{
Korkmaz M. C. and Chesneau, C., (2021). On the unit Burr-XII distribution with the quantile regression modeling and applications. \emph{Computational and Applied Mathematics}, \bold{40}(29), 1--26.
}
\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}

André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
