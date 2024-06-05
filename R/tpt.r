#-----------------------------------------------------------------------------#
#                                                                             #
#  MATRIX-BASED FLEXIBLE PROJECT PLANNING                                     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan, Aamir Saghir                                #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: June 2024                                                    #
#-----------------------------------------------------------------------------#

#' @export
tpt<- function(DSM,TD,SST=NULL)
{
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
      call. = FALSE
    )
  }
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

  DSM<-round(DSM)
  #T<-Re(T)
  N<-pracma::numel(TD)
  EST<-matrix(rep(0,N))         # EST will be an N by 1 null vector at the first step
  TD[matrix(diag(DSM)==0)*1]<-0 # The excluded task duration is irrelevant=>set to be 0
  TD<-matrix(TD)                # TD must  be column vector
  EFT<-EST+TD                   # EFTi=ESTi+Ti (i=1..N)
  DSM[(diag(DSM)==0)*(1:N),]<-0
  DSM[,(diag(DSM)==0)*(1:N)]<-0

  for (i in c(1:(N-1))){              # Forward pass
    for (j in c((i+1):N)){
      if (DSM[i,i]>0) {
        if (DSM[j,j]>0) {
          if (DSM[i,j]>0){
            if (EST[j]<EFT[i]){
              EST[j]<-EFT[i]
              EFT[j]<-EST[j]+TD[j]
            }
          }
        }
      }
    }
  }
  TPT<- max(EFT)                 # TPT is the makespan of the longest path
  LFT<- pracma::repmat(TPT,N,1)  # LFTi=TPT (i=1..N)
  LST<- LFT-TD                   # LSTi=LFTi-TDi (i=1..N)
  for (i in (N:2)){              # Backward pass
    for (j in ((i-1):1)){
      if (DSM[i,i]>0) {
        if (DSM[j,j]>0) {
          if (DSM[j,i]>0){
            if (LST[i]<LFT[j]){
              LFT[j]<-LST[i]
              LST[j]<-LFT[j]-TD[j]
            }
          }
        }
      }
    }
  }
  colnames(EST)<-"EST"
  colnames(EFT)<-"EFT"
  colnames(LST)<-"LST"
  colnames(LFT)<-"LFT"
  for (i in 1:N)
    if (LST[i]<EST[i])
      LST[i]<-EST[i]
  if (!is.null(rownames(DSM))){
    rownames(EST)<-rownames(EFT)<-rownames(LST)<-rownames(LFT)<-rownames(DSM)
  }else{
    rownames(EST)<-rownames(EFT)<-rownames(LST)<-rownames(LFT)<-rownames(DSM)<-paste("a",1:nrow(DSM),sep="_")
  }
  if (is.null(SST)){
    output <- list(TPT=TPT,EST=EST,EFT=EFT,LST=LST,LFT=LFT,SST=EST,SFT=EFT)
  }else{
    SST<-matrix(SST)
    ## CONT
    N<-pracma::numel(TD)
    TPT=0
    if (pracma::numel(which(diag(DSM)>0))){
      #Initialisation
      SFT=SST+TD
      #Forward pass
      for (i in c(1:(N-1))){
        for (j in c((i+1):N)){
          if (DSM[i,i]>0){
            if (DSM[j,j]>0){
              if (DSM[i,j]>0){
                if (SST[j]< SFT[i]){
                  SST[j]<-SFT[i]
                  SFT[j]<-SST[j]+TD[j]
                }
              }
            }
          }
        }
      }
      TPT<-max(SFT)
      colnames(SST)<-"SST"
      colnames(SFT)<-"SFT"
      if (!is.null(rownames(DSM))){
        rownames(SST)<-rownames(SFT)<-rownames(DSM)
      }else{
        rownames(SST)<-rownames(SFT)<-rownames(DSM)<-paste("a",1:nrow(DSM),sep="_")
      }
    }
    output <- list(TPT=TPT,EST=EST,EFT=EFT,LST=LST,LFT=LFT,SST=SST,SFT=SFT)
  }
  class(output)<-"TPT"
  return(output)
}
