\name{su}
\alias{su}
\title{
 Summary
}
\description{
 \code{su} produces summary for given percentiles
}
\usage{
su(x,
  mat.var.in.col = TRUE,
  digits = 4,
  probs = c(0.1, 0.25, 0.5, 0.75, 0.9),
  print = FALSE)
}
\arguments{
 \item{x}{
  an object 'data.frame', 'matrix' or 'vector'.
 }
 \item{mat.var.in.col}{
logical. If \code{TRUE} summary is applied by columns.
}
 \item{digits}{
number of digits.
}
\item{probs}{
vector of percentiles.
}
\item{print}{
 logical. If \code{TRUE}, summary is outputed on screen.
 }
}

\details{
 No details

}

\value{
 \code{su} returns a vector or matrix with summary
}

\author{
 Anonymouso@authors.com
}

\seealso{
\code{\link[base]{summary}}
}

\examples{
t2 <- as.matrix(rnorm(10))
t3 <- as.matrix(rnorm(10))
t2[2] <- NA
su(t2)
colnames(t2) <- "try1"
su(t2)
su(t(t2))
summary(t2)

su(rnorm(11))

my.matrix <- cbind(t2,t3)
colnames(my.matrix) <- c("t22", "t32")

summary(my.matrix)
su(my.matrix, T)

summary(data.frame(t21 = t2, t31 = t3))
su(data.frame(t21 = t2, t31 = t3))

}


\keyword{ Stochastic Frontier Analysis }
\keyword{ Heteroskedasticity }
\keyword{ Parametric efficiency analysis }

