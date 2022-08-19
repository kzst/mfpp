\name{tpq}
\alias{tpq}
\title{ Function to calcualte Total Project Quality for a project structure.
}
\description{ Calculate Total Project Quality of a project structure.}
\usage{
tpq (DSM,PEM,q,QD=NULL)
}

\arguments{
  \item{DSM}{N by N Upper triangular binary matrix of logic domain (a project structure
of a PEM matrix)(a numeric matrix).}
  \item{PEM}{N by N Upper triangular binary matrix of logic domain (a numeric matrix).}
  \item{q}{N by 1 vector of quality parameters (a numeric vector).}
  \item{QD}{N by w Quality Domain to calculate relative TPQ. The default value is QD=NULL}
}

\value{
 \item{TPQ}{ Total Project Quality (a scalar).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan, Aamir Saghir*

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpr}}, \code{\link{tpt}}, \code{\link{maxscore_PEM}}.}
\examples{

# Calculatation of total project quality of a project using MFPP package.

library(pracma)

# Generation of a 10 by 10 upper triangular binary matrix (PEM) of logic domain of a project.

PEM <- triu(rand(10)*0.5+0.5)

# Generation of a 10 by 10 Upper triangular binary matrix of logic domain
# (a project structureof a PEM matrix).

DSM <- round(PEM)
# Generation of a 10 by 1 vector of quality parameters of a project.
q <- rand(10,1)

# For fixed project plan the TPQ is the geometric mean of quality parameters

tpq(DSM,DSM,q)

# For flexibility project plan the TPQ
# is the weighted geometric mean of
# quality parameters

tpq(DSM,PEM, q)

QD2<-cbind(q,rand(10,1)) # Generate two completion modes

# Calculate relative quality considering the best choice

tpq(DSM,PEM,q,QD2)

}
\concept{project planning}
