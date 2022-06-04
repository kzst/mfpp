\name{tpt}
\alias{tpt}
\title{ Function to evaluate EST, EFT, LST and LFT times of activity of a project.
%%  ~~function to do ... ~~
}
\description{ Calculate EST, EFT, LST and LFT times of activity of a project.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tpt(DSM,TD,SST=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DSM}{N by N Upper triangular binary matrix of logic domain (a numeric matrix).}
  \item{TD}{N by 1 vector of time duration (a numeric vector).}
  \item{SST}{N by 1 vector of scheduled start time (a numeric vector).}
}

\value{
 \item{TPT}{ Total Project Time (a scalar).}
 \item{EST}{ Early Start Time (a vector).}
 \item{EFT}{ Early Finish Time (a vector).}
 \item{LST}{ Latest Start Time (a vector).}
 \item{LFT}{ Latest Finish Time (a vector).}
 \item{SST}{ Scheduled Start Time (a vector).}
 \item{SFT}{ Scheduled Finish Time (a vector).}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpr}}, \code{\link{maxscore_PEM}}.}
\examples{

# Calculatation of EST, EFT, LST and LFT times of activity of a project

# Define a 3 by 3 upper triangular binary matrix (DSM) of logic domain (LD) of a project.

# Specify sequence of three tasks
# Diagonal values indicates the priority of tasks
# 1=the task is mandatory
LD<-rbind(c(1,1,0), c(0,1,1), c(0,0,1))

# Define a 3 by 1 vector of time durations of a project, where the durations are 3,4,5 respectively.
TD<-c(3,4,5)

# Calculate project duration (total project time, TPT)
TPT<-tpt(LD,TD)

summary(TPT)

# Specify initial scheduled start time, where the delay is 1 for the .
SST <- c(1,0,0)

# Including scheduled starts and finishes.

TPT<-tpt(LD,TD,SST)

# Print schedule table
summary(TPT)

# Plot the scheduled Gantt chart
plot(TPT,sched="S")
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
