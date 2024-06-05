\name{percent}
\alias{percent}
\title{ Function to calculate desired project completion characteristic of a project structure.}
\description{ Function to calculate project completion desired charcateristic type either cost, or quality and/or resources with or withoud quality parameters, or score ratio or time duration of a project scenario for the given ratio and completion modes.}
\usage{
percent(PDM,type=c("c","q","qd","r","s","t"),w=2,Rs=2,ratio=1)
}

\arguments{
  \item{PDM}{is a Project Domain Matrix (a numeric matrix), where PDM=[LD,TD, CD, <QD,RD>]. LD is an N by N upper triangular matrix of logic domain (a numeric matrix); TD is an N by w matrix of task durations (a numeric matrix); CD is an N by w matrix of cost demands (a numeric matrix); QD is an optinal N by w matrix of quality parameters (a numeric optional matrix); RD is an N by w*nR matrix of resource demands (a numeric optional matrix).}
  \item{type}{ Type to calculate the project completion characteristic either cost or quality with or without defined quality parameters, or resources with or without quality parameters, or minimum completion score or minimum task duration of a project scenario (a logical value).}
   \item{w}{ Number of completion modes (a numeric integer).}
   \item{Rs}{ Number of resources (a numeric integer).}
   \item{ratio}{ratio of desired characteristic of a project between interval [0,1] (a numeric value).}
}

\value{
 \item{w}{ Number of completion modes (a numeric integer).}
 \item{Rs}{ Number of resources (a numeric integer).}
 \item{ratio}{ratio of desired characteristic of a project between interval [0,1] (a numeric value).}
 \item{Cc}{Project completion cost (a numeric value).}
 \item{Cq}{Project completion quality when quality parameters are not considered (a numeric value).}
 \item{Cq}{Project completion quality when quality parameters are considered (a numeric value).}
 \item{CR}{Minimum resources for the completion when quality parameters are not considered (a numeric value).}
  \item{CR}{Minimum resources for the completion when quality parameters are considered (a numeric value).}
   \item{Cs}{Minimum completion score ratio of a project (a numeric value).}
  \item{Ct}{Minimum task duration to project a project (a numeric value).}
}
\references{
Kosztyán, Z. T. (2022). MFPP: Matrix-based flexible project planning. SoftwareX, 17, 100973.
}
\author{ Zsolt T. Kosztyan*, Aamir Saghir

e-mail: kzst@gtk.uni-pannon.hu}

\seealso{\code{\link{tpc}}, \code{\link{tpq}}, \code{\link{tpt}}, \code{\link{tpr}},\code{\link{maxscore_PEM}}, \code{\link{generatepdm}}.}
\examples{
## Examples for verification
w<-2 # Number of completion modes
Rs<-2 # Number of resources
ratio<-1.0 # Percentage 0-1

# Specification of Logic Domain
LD<-matrix(c(1,1,0,1,0,0,
             0,0.8415,0,0,0,0,
             0,0,1,1,0,0,
             0,0,0,0.6478,0,0,
             0,0,0,0,0,0,
             0,0,0,0,0,0),nrow=6, byrow=TRUE)

# Specification of Time Domain
TD<-matrix(c(1.9755,2.3408,
             2.6564,2.7002,
             0.9018,1.1077,
             1.0237,1.2117,
             0,0,
             0,0),nrow=6, byrow=TRUE)
# Specification of Cost Domain
CD<-matrix(c(2.0768,2.3008,
             0.4819,0.5279,
             3.4315,3.7682,
             3.5902,3.8245,
             0,0,
             0,0),nrow=6, byrow=TRUE)

# Specification of Cost Domain
QD<-matrix(c(0.8723,0.9102,
             0.9017,0.8215,
             0.7250,0.9911,
             1.0000,0.8245,
             0,0,
             0,0),nrow=6, byrow=TRUE)


# Specification of Resource Domain
RD<-matrix(c(2.2130,2.4636,1.7948,2.0603,
             1.9709,2.1952,0.8216,0.8890,
             1.6255,1.9432,1.8715,2.2341,
             1.1377,1.3528,1.6768,1.8769,
             0,0,0,0,
             0,0,0,0),nrow=6, byrow=TRUE)

# PDM without QD
PDM<-cbind(LD,TD,CD,RD)

# Percentiles without QD
CONST<-percent(PDM,type=c("c","r","s","t"),w,Rs,ratio)
CONST

# PDM with QD
PDM<-cbind(LD,TD,CD,QD,RD)

# Percentiles with QD
CONST<-percent(PDM,type=c("c","q","r","s","t"),w,Rs,ratio)
CONST

}
\concept{project planning}
