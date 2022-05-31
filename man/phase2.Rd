\name{phase2}
\alias{phase2}
\title{ Function to simulate shock effects.
%%  ~~function to do ... ~~
}
\description{ Simulating the shock effects of the changes of customer claims. Two steps of Monte-Carlo simulation is applied for selecting (1) and changing (2) task demands and the structures.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phase2 (x,p=0.10, s=5.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{is either a Project Domain Matrix (a numeric matrix), where PDM=[LD,TD, CD, {QD,RD}]. LD is an N by N upper triangular matrix of logic domain (a numeric matrix); TD is an N by w matrix of task durations (a numeric matrix); CD is an N by w matrix of cost demands (a numeric matrix); QD is an optinal N by w matrix of quality parameters (a numeric optional matrix); RD is an N by w*nR matrix of resource demands (a numeric optional matrix), or PDM list, which contain a PDM matrix, and the number of completion modes (w) and the number of resources.}
   \item{p}{ Probability factor for task selction (a numeric value between 0 and 1). The default value is p=1.00}
   \item{s}{ Scale factor: the ratio of the modification (a numeric integer). The default value is s=5.0}
}

\value{
 \item{PDMout}{PDM matrix with same structure as the input PDM matrix (a matrix).}

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

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{tpr}},\code{\link{maxscore_PEM}},\code{\link{phase1}},\code{\link{phase3}}, \code{\link{generatepdm}}.}
\examples{

# Simulation of project domain matrix based on customer claims
# of a project scenario using MFPP package.

# Generate a project domain matrix based on
# logic domain, task duration, cost demans,
# quality parameters (optional),  resource demand (optional) of a project.

PDM=generatepdm(30,0.05,0,20,30,20,2,2,2,QD=TRUE,lst=TRUE)
PDM

# Define probability factor for task selection of a project.
p=0.10

# Choose a scale factor for the ration of modification.
s=5.0

# Simulate the shock effect of a project scenario using MFPP package
# with default values.
 Res<-phase2(PDM)
 summary(Res)

# Simulate the effect of a project scenario using MFPP package
# with specified values.
 Res<-phase2(PDM$PDM,p,s)
 summary(Res,w=2,Rs=2)
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
