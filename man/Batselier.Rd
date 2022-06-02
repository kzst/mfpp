\name{Batselier}
\alias{Batselier}
\docType{data}
\title{
Real-life project database by Batselier and Vanhoucke (2015)
}
\description{
Real-life project database of single-mode individual projects by Batselier and Vanhoucke (2015)
}
\usage{data("Batselier")}
\format{
  A list (Collection_PDM) of PDM_list and PDM_constraint types of lists, where
  \describe{
    \item{\code{PDM_list}}{contains PDM_matrix (PDM)  and the number of completion modes (w), and number of resources (Rs)}
    \item{\code{PDM_const}}{a list of (renewable resource) constraints}
  }
}

\references{
 Batselier, J. & Vanhoucke, M. (2015). Construction and evaluation framework for a real-life project database. International Journal of Project Management, 33(3), 697–710, https://doi.org/10.1016/j.ijproman.
}
\examples{
data(Batselier)

# Summary of data collection
summary(Batselier)

# Summary of the first project structure
summary(Batselier[[1]]$PDM_list)

# Summary of the first project constraints
summary(Batselier[[1]]$PDM_const)
}
\keyword{datasets}
