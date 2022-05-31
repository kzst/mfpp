\name{tpq}
\alias{tpq}
\title{ Function to calcualte Total Project Quality for a project structure.
%%  ~~function to do ... ~~
}
\description{ Calculate Total Project Quality of a project structure.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tpq (DSM,PEM,q,QD=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DSM}{N by N Upper triangular binary matrix of logic domain (a project structure
of a PEM matrix)(a numeric matrix).}
  \item{PEM}{N by N Upper triangular binary matrix of logic domain (a numeric matrix).}
  \item{q}{N by 1 vector of quality parameters (a numeric vector).}
  \item{QD}{N by w Quality Domain to calculate relative TPQ. The default value is QD=NULL}
}

\value{
 \item{TPQ}{ Total Project Quality (a scalar).}

%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan, Aamir Saghir*

e-mail: kzst@gtk.uni-pannon.hu
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

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

# Calculation of total project quality using MFPP package.
tpq(DSM,PEM, q)


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
