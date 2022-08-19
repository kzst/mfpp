\name{maxscore_PEM}
\alias{maxscore_PEM}
\title{ Function to calculate maximal score value (PMAX) of possible project scenarios.
}
\description{ Calculate maximal score value (PMAX) of possible project scenarios.
}
\usage{
maxscore_PEM(PEM, P=PEM, Q=1-PEM)
}
\arguments{
  \item{PEM}{N by N upper triangular adjacency matrix of logic network (a numeric matrix).}
  \item{P}{N by N score matrix of task/dependency inclusion (a numeric matrix). The default P matrix is P=PEM}
  \item{Q}{N by N score matrix of task/dependency exclusion (a numeric matrix). The default Q matrix is Q=1-PEM}
}

\value{
 \item{score}{ The maximal score value of the project scenario (a scalar).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpr}}, \code{\link{tpt}}.}
\examples{

# Calculatation of the maximal score value of the project scenario using MFPP package.


# Define a 3 by 3 upper triangular adjacency matrix (PEM) of logic domain of a project.

PEM <- rbind(c(0.8,0.4,0.8),
             c(0.0,0.7,0.7),
             c(0.0,0.0,0.4))

# Define a 3 by 3 score matrix of task/dependency inclusion.
P <- PEM

# Define a 3 by 3 score matrix of task/dependency exclusion.
Q <- 1-P

# Calculation of the maximal score value of the project using MFPP package.
maxscore_PEM(PEM,P, Q)
}
\concept{project planning}
