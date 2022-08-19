\name{tpr}
\alias{tpr}
\title{ Function to calculate maximum resource demands of a project.}
\description{ Calculate maximum resource demands of a project scenario.}
\usage{
tpr (SST,DSM,TD, RD,res.graph=FALSE)
}
\arguments{
 \item{SST}{N by 1 vector of Scheduled Start Time (a numeric vector).}
  \item{DSM}{N by N binary upper triangular matrix of the logic domain (a numeric matrix).}
  \item{TD}{ N by 1 vector of task durations (a numeric vector).}
  \item{RD}{ N by nR matrix of resource demands (a numeric matrix).}
  \item{res.graph}{Logic value of whether plot or not the resource graph}
}
\value{
 \item{rMAX}{an nR by 1 vector of maximum resource demands (a scalar vector).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{maxscore_PEM}}.}
\examples{

# Calculatation of maximum resource demands of a project scenario using MFPP package.

library(pracma)

# Generate a 10 by 10 upper triangular binary matrix (DSM) of logic domain of a project.

DSM <- triu(round(rand(10)))

# Generate a 10 by 1 vector of task durations of a project.
TD <- rand(10,1)*20

# Define a 10 by 1 vector of Scheduled Start Time using tpt function.
SST <- tpt(DSM,TD)[["EST"]]

# Generate a 10 by 3 matrix of resource demands of a project.
RD <- rand(10,3)*5

# Calculation of maximum resource demands of a project using MFPP package.
 tpr(SST,DSM,TD, RD)

## Plot resources for SST
tpr(SST,DSM,TD,RD,res.graph = TRUE)
}
\concept{project planning}
