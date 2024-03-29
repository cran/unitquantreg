% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpqr-leeg.R
\name{leeg}
\alias{leeg}
\alias{dleeg}
\alias{pleeg}
\alias{qleeg}
\alias{rleeg}
\title{The Log-extended exponential-geometric distribution}
\usage{
dleeg(x, mu, theta, tau = 0.5, log = FALSE)

pleeg(q, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

qleeg(p, mu, theta, tau = 0.5, lower.tail = TRUE, log.p = FALSE)

rleeg(n, mu, theta, tau = 0.5)
}
\arguments{
\item{x, q}{vector of positive quantiles.}

\item{mu}{location parameter indicating the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.}

\item{theta}{nonnegative shape parameter.}

\item{tau}{the parameter to specify which quantile is to be used.}

\item{log, log.p}{logical; If TRUE, probabilities p are given as log(p).}

\item{lower.tail}{logical; If TRUE, (default), \eqn{P(X \leq{x})} are returned, otherwise \eqn{P(X > x)}.}

\item{p}{vector of probabilities.}

\item{n}{number of observations. If \code{length(n) > 1}, the length is taken to be the number required.}
}
\value{
\code{dleeg} gives the density, \code{pleeg} gives the distribution function,
\code{qleeg} gives the quantile function and \code{rleeg} generates random deviates.

Invalid arguments will return an error message.
}
\description{
Density function, distribution function, quantile function and random number generation function
for the Log-extended exponential-geometric distribution reparametrized in terms of the \eqn{\tau}-th quantile, \eqn{\tau \in (0, 1)}.
}
\details{
Probability density function
\deqn{f(y\mid \alpha ,\theta )=\frac{\theta \left( 1+\alpha \right) y^{\theta -1}}{\left( 1+\alpha y^{\theta }\right) ^{2}}}

Cumulative distribution function
\deqn{F(y\mid \alpha ,\theta )=\frac{\left( 1+\alpha \right) y^{\theta }}{1+\alpha y^{\theta }}}

Quantile function
\deqn{Q(\tau \mid \alpha ,\theta )=\left[ \frac{\tau }{1+\alpha \left( 1-\tau\right) }\right] ^{\frac{1}{\theta }}}

Reparameterization
\deqn{\alpha=g^{-1}(\mu )=-\frac{1-\tau \mu ^{\theta }}{\left( 1-\tau \right) }}
}
\examples{
set.seed(123)
x <- rleeg(n = 1000, mu = 0.5, theta = 1.5, tau = 0.5)
R <- range(x)
S <- seq(from = R[1], to = R[2], by =  0.01)
hist(x, prob = TRUE, main = 'Log-extended exponential-geometric')
lines(S, dleeg(x = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(ecdf(x))
lines(S, pleeg(q = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
plot(quantile(x, probs = S), type = "l")
lines(qleeg(p = S, mu = 0.5, theta = 1.5, tau = 0.5), col = 2)
}
\references{
Jodrá, P. and Jiménez-Gamero, M. D., (2020). A quantile regression model for bounded responses based on the exponential-geometric distribution. \emph{Revstat - Statistical Journal}, \bold{18}(4), 415--436.
}
\author{
Josmar Mazucheli \email{jmazucheli@gmail.com}

André F. B. Menezes \email{andrefelipemaringa@gmail.com}
}
