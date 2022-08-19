\name{get.structures}
\alias{get.structures}
\title{ Function to calculate minimal/maximal/most likely project structures.
}
\description{ Specify minimal, maximal, maximin, minimax, and most likely/most desired project structures.
}
\usage{
get.structures(x,type=c("min","max","minimax","maximin","most"))
}
\arguments{
  \item{x}{N by M PDM matrix or a PDM list.}
  \item{type}{type of specified project structure. "min" is a minimal strucuture, where the time, cost demands, and as well as the quality parameters are minimal; "max" where the time, cost demands, and as well as the quality parameters are maximal; "minimax", where the cost and resource demands are minimal; "maximix", where the cost and resource demands are maxnimal; "most", where the total project score is maximal.}
}

\value{
 \item{PDMs}{set of PDM_matrix or a set of PDM_list (depending on the input type)}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}

\seealso{\code{\link{generatepdm}}.}
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
PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4)

# Get main structures.

PDMs<-get.structures(PDM)
PDMs

# Generation of PDM list for TCTP a project scenario using MFPP package.
PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=TRUE)

# Get main structures.

PDMs<-get.structures(PDM)
PDMs

}
\concept{project planning}
