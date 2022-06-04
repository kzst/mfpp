\name{tpc}
\alias{tpc}
\title{ Function of Cost demands of a project
%%  ~~function to do ... ~~
}
\description{ Calculate cost demands of a project structure.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tpc(DSM,CD)
}
\arguments{
  \item{DSM}{N by N Upper triangular binary matrix of logic domain (a numeric matrix).}
  \item{CD}{N by 1 vector of cost demands (a numeric vector).}
}

\value{
 \item{TPC}{ Total Project Cost (a scalar).}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
%%  ~~who you are~~
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{tpc}}, \code{\link{tpt}}, \code{\link{tpq}}, \code{\link{tpr}}, \code{\link{maxscore_PEM}}.}
\examples{

# Calculation of TPC for a generated project scenario.
# using MFPP package.

# Generation of 10 by 10 Upper triangular binary matrix (DSM) of logic domain.

library(Matrix)
library(pracma)

DSM<-triu(round(rand(10)*.5+.5)) # Generate DSM

# Generation of 10 by 1 cost vector (cost demands 30)
C<-rand(10,1)*30                 # Generate C vector (cost demands)


# Calculation of total project cost using MFPP package.
TPC<-tpc(DSM,C)                  # Calculate TPC



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
