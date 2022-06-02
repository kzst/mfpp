\name{Boctor}
\alias{Boctor}
\docType{data}
\title{
Sumulated project database by Boctor (1993)
}
\description{
Sumulated project database of multimode, individual projects by Boctor (1993)
}
\usage{data("Boctor")}
\format{
  A list (Collection_PDM) of PDM_list and PDM_constraint types of lists, where
  \describe{
    \item{\code{PDM_list}}{contains PDM_matrix (PDM)  and the number of completion modes (w), and number of resources (Rs)}
    \item{\code{PDM_const}}{a list of (renewable resource) constraints}
  }
}

\references{
Boctor, F. F. (1993). Heuristics for scheduling projects with resource restrictions and several resource-duration modes. The international journal of production research, 31(11), 2547–2558, https://doi.org/10.1080/00207549308956882.
}
\examples{
data(Boctor)

# Summary of data collection
summary(Boctor)

# Summary of the first project structure
summary(Boctor[[1]]$PDM_list)

# Summary of the first project constraints
summary(Boctor[[1]]$PDM_const)
}
\keyword{datasets}
