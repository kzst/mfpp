\name{minscore_PEM}
\alias{minscore_PEM}
\title{ Function to calculate minimal score value of possible project scenarios.}
\description{ Calculate minimal score value (PMIN) of possible project scenarios.}
\usage{
minscore_PEM(PEM,P=PEM, Q=1-PEM)
}
\arguments{
  \item{PEM}{An N by N upper triangular adjacency matrix of logic network (a numeric matrix).}
  \item{P}{ N by N score matrix of task/dependency inclusion (a numeric matrix). The default P matrix is P=PEM}
  \item{Q}{ N by N score matrix of task/dependency exclusion (a numeric matrix). The default Q matrix is Q=1-PEM}
}

\value{
 \item{score}{minimal score value of a project (a numeric value).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{tpr}},\code{\link{maxscore_PEM}}, \code{\link{generatepdm}}, \code{\link{phase3}}, \code{\link{percent}}.}
\examples{

# Calculation of minimal score value of a project scenarios using MFPP package.

# Define a 3 by 3 upper traingular logic domain of a project scenario.

PEM= matrix(c(0.8,0.4,0.8,
              0.0,0.7,0.7,
              0.0,0.0,0.4), ncol=3, byrow=TRUE)

# Define a score matrix of task/dependency inclusion.
P=PEM
# Define a score matrix of task/dependency exclusion.
Q=1-P

# Calculation of minimal score value for the above specification
# of a project scenario using MFPP package.

minscore_PEM(PEM,P,Q)
}
\concept{project planning}
