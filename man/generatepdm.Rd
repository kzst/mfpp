\name{generatepdm}
\alias{generatepdm}
\title{ Function to generate a PDM matrix for flexible project planning problems.
}
\description{ Function to generate a PDM matrix for flexible project planning problems.
}
\usage{
generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=FALSE)
}

\arguments{
  \item{N}{ Number of tasks (a numeric value)}
  \item{ff}{flexibility factor between [0,1] (a numeric value).}
  \item{cf}{connectivity factor (a numeric integer).}
  \item{mTD}{max value of TD (a numeric positive integer).}
  \item{mCD}{max value of CD (a numeric positive integer).}
  \item{mRD}{max value of RD (a numeric positive integer).}
  \item{w}{number of modes  (a numeric integer).}
  \item{nR}{number of resources (a numeric integer).}
   \item{nW}{number of possible extra tasks (a positive numeric integer).}
    \item{scale}{scale parameter, the default value is 1.4 (a numeric value).}
   \item{QD}{logic parameter, the default value is FALSE, which means, the quality domain is neglected.}
   \item{lst}{(logic parameter, the default value is FALSE, which means, the output is a PDM matrix, while TRUE means the output is a PDM list including number of modes (w) and the number of resources (Rs).}
}

\value{
  \item{PDM}{[PEM,TD, CD, <QD,RD>] is a Project Domain Matrix (a numeric matrix), where;}
  \item{PEM}{An N by N upper triangular matrix of logic domain (a numeric matrix).}
  \item{TD}{N by w matrix of task durations (a numeric matrix).}
  \item{CD}{N by w matrix of cost demands (a numeric matrix).}
  \item{QD}{N by w matrix of quality parameters (a numeric matrix, optional).}
  \item{RD}{ N by w*nR matrix of resource demands (a numeric matrix, optional).}
  \item{w}{Number of completion modes (scalar, if the output is PDM list (lst=TRUE)).}
  \item{Rs}{ Number of resources (scalar, if the output is PDM list (lst=TRUE)).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{tpr}},\code{\link{maxscore_PEM}}, \code{\link{generatepdm}}.}
\examples{

# Generation of PDM matrix for flexible project planning MFPP package.

# Define number of modes, flexibility factor and connectivity factor of a project scenerio.
N=4;ff=0.05;cf=0

# Define maximum value of time domain, Cost domain and Resourcces domain of a project scenerio.
mTD=3;mCD=4;mRD=3

# Define number of modes, number of resources,
# number of possible extra tasks, scale and quality domain of a project scenerio.

w=2;nR=2;nW=2
scale=1.4  #(default value)

# Generation of PDM matrix for TCTP a project scenario using MFPP package.
generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4)

# Generation of PDM list for TCTP a project scenario using MFPP package.
generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=TRUE)


# Generation of PDM matrix for TQCTP a project scenario using MFPP package.
generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=TRUE)

 # Generation of PDM list for TQCTP a project scenario using MFPP package.
generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=TRUE,lst=TRUE)

}

\concept{project planning}

