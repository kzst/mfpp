\name{tpc}
\alias{tpc}
\title{ Function of Cost demands of a project.}
\description{ Calculate cost demands of a project structure.}
\usage{
tpc(DSM,CD)
}
\arguments{
  \item{DSM}{N by N Upper triangular binary matrix of logic domain (a numeric matrix).}
  \item{CD}{N by 1 vector of cost demands (a numeric vector).}
}

\value{
 \item{TPC}{ Total Project Cost (a scalar).}}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}
\seealso{\code{\link{tpc}}, \code{\link{tpt}}, \code{\link{tpq}}, \code{\link{tpr}}, \code{\link{maxscore_PEM}}.}
\examples{

# Calculation of TPC for a generated project scenario.
# using MFPP package.

# Generation of 10 by 10 Upper triangular binary matrix (DSM) of logic domain.

library(Matrix)
library(pracma)

DSM<-triu(round(rand(10)*.5+.5)) # Generate DSM

# Generation of 10 by 1 cost vector (cost demands 30)
C<-rand(10,1)*30                 # Generate C vector (cost demands)


# Calculation of total project cost using MFPP package.
TPC<-tpc(DSM,C)                  # Calculate TPC

}

\concept{project planning}
