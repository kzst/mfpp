\name{phase3}
\alias{phase3}
\title{ Function to simulate the effects of the change of customer claims.}
\description{ Phase 3 selects P percent of nodes (i.e. tasks) or arcs (i.e. dependencies), where the score value can be changed up to the maximal change effects.}
\usage{
phase3 (x,p=0.10, s=0.50, nW=0)
}
\arguments{
  \item{x}{is either a Project Domain Matrix (a numeric matrix), where PDM=[LD,TD, CD, <QD,RD>]. LD is an N by N upper triangular matrix of logic domain (a numeric matrix); TD is an N by w matrix of task durations (a numeric matrix); CD is an N by w matrix of cost demands (a numeric matrix); QD is an optinal N by w matrix of quality parameters (a numeric optional matrix); RD is an N by w*nR matrix of resource demands (a numeric optional matrix), or PDM list, which contain a PDM matrix, and the number of completion modes (w) and the number of resources.}
   \item{p}{ Probability of structural change (default value is p=0.10).}
   \item{s}{ Increase/decrease ratio of priorities (default value is s=0.50).}
   \item{nW}{ Number of unplanned tasks (default value is nW=0).}
}

\value{
 \item{PDMout}{PDM matrix with same structure as the input PDM matrix (a matrix).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{tpr}},\code{\link{phase1}},\code{\link{phase2}},\code{\link{maxscore_PEM}}, \code{\link{generatepdm}}.}
\examples{

# Simulation of project domain matrix based on customer claims
# of a project scenario using MFPP package.

# Generate a project domain matrix based on
# logic domain, task duration, cost demans,
# quality parameters (optional),  resource demand (optional) of a project.

PDM=generatepdm(30,0.05,0,20,30,20,2,2,2,QD=TRUE,lst=TRUE)
PDM

# Define probability factor for task selection and dependency selecton of a project.
p=0.10

# Choose a scale factor for the ration of modification.
s=0.10

# Simulate the effects of change of customer claims of a project scenario using MFPP package
# with default values.

 Res<-phase3(PDM)
 summary(Res)

# Simulate the effects of change of customer claims of a project scenario using MFPP package
# with specified values.
 Res<-phase3(PDM$PDM,p,s,nW=2)
 summary(Res,w=2,Rs=2)
}
\concept{project planning}
