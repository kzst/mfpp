\name{plot.mfpp}
\alias{plot.PDM_matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot function for Matrix-Based Flexible Project Planning
}
\description{
Plot function for the Package of Matrix-Based Flexible Project Planning (mfpp)
}
\usage{
\method{plot}{PDM_matrix}(x,w=getOption("w"),Rs=getOption("Rs"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class 'PDM_matrix'.}
  \item{w}{number of completion modes (positive integer)}
  \item{Rs}{number of resources (positive integer)}
  \item{...}{other graphical parameters.}
}

\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{
Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}


\seealso{\code{\link{summary}}, \code{\link{generatepdm}}, \code{\link{percent}}.}

\examples{
# Plot PDM matrix.

# Define number of modes, flexibility factor and connectivity factor of a project scenerio.
N=4;ff=0.05;cf=0

# Define maximum value of time domain, Cost domain and Resourcces domain of a project scenerio.
mTD=3;mCD=4;mRD=3

# Define number of modes, number of resources,
# number of possible extra tasks, scale and quality domain of a project scenerio.

w=2;nR=2;nW=2
scale=1.4  #(default value)

# Generation of PDM list for TCTP a project scenario using MFPP package.

PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4)

plot(PDM)
}

\concept{project management}
