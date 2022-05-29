\name{summary.mfpp}
\alias{summary.PDM_const}
\alias{summary.PDM_matrix}
\alias{summary.PDM_list}
\alias{summary.Set_PDM_matrix}
\alias{summary.Set_PDM_list}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print PDM constraints, matrices, lists, sets.
}

\description{
Print PDM constraints, matrices, lists, sets.
}

\usage{
\method{summary}{PDM_const}(object, digits =  getOption("digits"), ...)
\method{summary}{PDM_matrix}(object, digits =  getOption("digits"),w=getOption("w"),Rs=getOption("Rs"), ...)
\method{summary}{PDM_list}(object, digits =  getOption("digits"), ...)
\method{summary}{Set_PDM_matrix}(object, digits =  getOption("digits"),w=getOption("w"),Rs=getOption("Rs"), ...)
\method{summary}{Set_PDM_list}(object, digits =  getOption("digits"), ...)

}

%- maybe also 'usage' for other xs documented here.
\arguments{
  \item{object}{an object of class 'PDM_const', 'PDM_matrix', 'PDM_list', 'Set_PDM_matrix', or 'Set_PDM_list'.}
  \item{digits}{the number of significant digits to use when \code{add.stats = TRUE}.}
  \item{w}{number of completion modes (positive integer)}
  \item{Rs}{number of resources (positive integer)}
  \item{...}{additional arguments affecting the print produced.}
}

\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{
Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}

\seealso{\code{\link{generatepdm}},\code{\link{get.structures}}, \code{\link{percent}},\code{\link{plot}}.}

\examples{

# Print PDM constriants, matrices, lists, and sets.

# Define number of modes, flexibility factor and connectivity factor of a project scenerio.
N=4;ff=0.05;cf=0

# Define maximum value of time domain, Cost domain and Resourcces domain of a project scenerio.
mTD=3;mCD=4;mRD=3

# Define number of modes, number of resources,
# number of possible extra tasks, scale and quality domain of a project scenerio.

w=2;nR=2;nW=2
scale=1.4  #(default value)

# Generation of PDM list for TCTP a project scenario using MFPP package.

PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,lst=TRUE)

# Generate PDM maximal constraints (ratio=1.0)

CONST<-percent(PDM)

# Calculate main structures

PDMs<-get.structures(PDM)

# print of PDM constraint

summary(CONST)

# print of PDM list

summary(PDM)

# print of main structures of PDM matrices

summary(PDMs)
}


\concept{project management}
