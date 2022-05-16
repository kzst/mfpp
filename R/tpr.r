#-----------------------------------------------------------------------------#
#                                                                             #
#  MATRIX-BASED FLEXIBLE PROJECT PLANNING                                     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan, Aamir Saghir                                #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: May 2022                                                     #
#-----------------------------------------------------------------------------#

#' @export
tpr<- function(SST,DSM,TD,RD)
{
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop(
      "Package \"pracma\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Rfast", quietly = TRUE)) {
    stop(
      "Package \"Rfast\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
      call. = FALSE
    )
  }
  N <-pracma::numel(SST)               #Number of tasks
  DSM <-round(pracma::triu(DSM))       #DSM must be an upper triangular binary matrix.
  DSM[(matrix(diag(DSM)==0)*1),]<-0    #Excluded task has no dependency
  DSM[,(matrix(diag(DSM)==0)*1)]<-0    #Excluded task has no dependency
  TD[matrix(diag(DSM)==0)*1]<-0        #Excluded task has no time demands
  RD[(matrix(diag(DSM)==0)*1),]<-0     #Excluded task has no resources
  RD[is.nan(RD)]<-0
  SST[matrix(diag(DSM)==0)*1]<-0       #Excluded task's start time is zero
  SST<-matrix(SST)                     #SST should be a column vector
  n<-pracma::numel(TD)                 #Number of elements in vector T
  SFT<-SST+TD                          #Initialisation
  for (i in 1: (n-1)) {
    for (j in (i+1):n) {
      if (DSM[i,i]> 0)
        if (DSM[j,j]> 0)
          if (DSM[i,j]> 0)             #If there is a dependency between task i and task j.
            if (SST[j] < SFT[i])
              SST[j]= SFT[i]
      SFT[j]= SST[j]+TD[j]
    }
  }
  BP<- sort(matrix(union(SST,SFT)))    #Breakpoints, where the resource demands should be recalculated
  b <- pracma::numel(BP)               # Number of breakpoints
  B<- t((pracma::repmat(SST,1,b)<=
           pracma::repmat(BP, N,1))& (pracma::repmat(SFT,1,b)>
                                        pracma::repmat(BP, N,1)))*1
  RESFUNC<- matrix(as.double(B), dim(B)[1], dim(B)[2])%*% RD       #RESFUNC=mtimes(B,R); %Calculate resource function
  rMAX<- Rfast::colMaxs(RESFUNC, value=TRUE)
  H<-t(rMAX)
  colnames(H)<-paste("R",1:ncol(H),sep="_")
  rownames(H)<-"TPR"
  return(H)
}
