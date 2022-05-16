\name{paretores}
\alias{paretores}
\title{ Calculate Pareto-optimal resource allocation
%%  ~~function to do ... ~~
}
\description{ Calculate Pareto-optimal resource allocation
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
paretores(DSM,TD,RD)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DSM}{An N by N upper triangular adjacency matrix of logic network with fixed dependencies (a binary matrix).}
  \item{TD}{ N by 1 matrix of task durations (a numeric matrix)}
  \item{RD}{ N by r matrix of task resources, where the number of resources is r.}
}

\value{
 \item{RD}{Pareto-optimal TPR (1 by r numeric vector)).}
 \item{SST}{Scheduled Start Time (N by N numeric matrix)).}

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

\seealso{\code{\link{tpr}}, \code{\link{percent}}.}
\examples{

# Calculation of Pareto-optimal resource allocation using MFPP package.

# Define a 3 by 3 upper traingular logic domain of a project structure.

# Specification of Logic Domain
DSM<-matrix(c(1,1,0,1,0,0,
              0,1,0,0,0,0,
              0,0,1,1,0,0,
              0,0,0,1,0,0,
              0,0,0,0,0,0,
              0,0,0,0,0,0),nrow=6, byrow=TRUE)

# Specification of Time Domain
TD<-matrix(c(2,
             3,
             1,
             4,
             0,
             0),nrow=6, byrow=TRUE)


# Specification of Resource Domain
RD<-matrix(c(4,2,
             3,2,
             5,1,
             6,4,
             0,0,
             0,0),nrow=6, byrow=TRUE)

RES<-paretores(DSM,TD,RD)
RES
# Verify results
tpr(RES$SST,DSM,TD,RD)

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
