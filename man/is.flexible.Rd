\name{is.flexible}
\alias{is.flexible}
\title{ Function to check the  flexibility of PDM matrix.
%%  ~~function to do ... ~~
}
\description{ Checking flexibility of PDM matrix or a PDM list.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
is.flexible(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{N by M PDM matrix or a PDM list.}
}


\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

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

# Flexibility check.

is.flexible(PDM)

# Flexibility check of rounded values.

is.flexible(round(PDM))

# Generation of PDM list for TCTP a project scenario using MFPP package.
PDM<-generatepdm(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=TRUE)

# Flexibility check.

is.flexible(PDM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{project planning}
% \concept{ ~cpt2 }
% Use only one concept per line.
