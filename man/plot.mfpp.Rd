\name{plot.mfpp}
\alias{plot.PDM_matrix}
\alias{plot.PDM_list}
\alias{plot.Set_PDM_matrix}
\alias{plot.Set_PDM_list}
\alias{plot.TPT}

\title{
Plot function for Matrix-Based Flexible Project Planning
}
\description{
Plot function for the Package of Matrix-Based Flexible Project Planning (mfpp)
}
\usage{
\method{plot}{PDM_matrix}(x,w=NULL,Rs=NULL,
type=c("orig","max","min","maximin","minimax","most","const"),
main=NULL,col=NULL,...)
\method{plot}{PDM_list}(x,
type=c("orig","max","min","maximin","minimax","most","const"),
main=NULL,col=NULL, ...)
\method{plot}{Set_PDM_matrix}(x,w=NULL,Rs=NULL,
type=c("orig","max","min","maximin","minimax","most","const"),col=NULL, ...)
\method{plot}{Set_PDM_list}(x,
type=c("orig","max","min","maximin","minimax","most","const"),col=NULL, ...)
\method{plot}{TPT}(x,
sched="E", ...)
}
\arguments{
  \item{x}{an object of class 'PDM_matrix' or 'PDM_list'.}
  \item{w}{number of completion modes (positive integer)}
  \item{Rs}{number of resources (positive integer)}
  \item{type}{set of plots, where (orig) means original, (min) minimal, (max) maximal, (maximin) maximin, (minimax) minimax structures are plotted, respectively, and as well as the minimal/maximal constraints (const)}
  \item{sched}{character starts are calculated either "E" to early, or "L" to late, or "S" to scheduled start time}
  \item{main}{title of the plot of fixed structure}
  \item{col}{color of barplot}
  \item{...}{other graphical parameters.}
}

\value{No return value, called for side effects}

\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{
Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}


\seealso{\code{\link{summary}}, \code{\link{generatepdm}}, \code{\link{percent}}.}

\examples{

# Plot PDM lists.

# Define number of modes, flexibility factor and connectivity factor of a project scenerio.
N=4;ff=0.3;cf=0

# Define maximum value of time domain, Cost domain and Resourcces domain of a project scenerio.
mTD=3;mCD=4;mRD=3

# Define number of modes, number of resources,
# number of possible extra tasks, scale and quality domain of a project scenerio.

w=2;nR=2;nW=2
scale=1.6

# Generation of PDM list for TCTP a project scenario using MFPP package.

PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale,lst=TRUE)

plot(PDM)
}

\concept{project management}
