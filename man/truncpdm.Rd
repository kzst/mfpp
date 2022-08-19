\name{truncpdm}
\alias{truncpdm}
\title{ Function to drop excluded tasks.}
\description{ Drop excluded tasks, where the diagonal value is 0.}
\usage{
truncpdm(x)
}
\arguments{
  \item{x}{N by M PDM matrix or a PDM list.}
}

\value{
 \item{PDM}{ N by M PDM matrix or a PDM list}
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

# Drop excluded tasks.

PDM<-truncpdm(PDM)

# Generation of PDM list for TCTP a project scenario using MFPP package.
PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=TRUE)

# Drop excluded tasks.

PDM<-truncpdm(PDM)

}
\concept{project planning}
