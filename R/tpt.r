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
tpt<- function(DSM,TD)
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
  dsm<-Matrix::triu(DSM,1)
  dsm[(matrix(diag(DSM)==0)*1),]<-0
  #dsm[,(matrix(diag(DSM)==0)*1)]<-0
  for (i in 1: N)          # Forward pass
  {
    est<-pracma::arrayfun("*",as.matrix(dsm),EFT)
    EST<-matrix(Rfast::colMaxs(est, value=TRUE)) #EST calculated step by step
    EFT[i]<-EST[i]+TD[i]
  }

  TPT<- max(EFT)       # TPT is the makespan of the longest path
  LFT<- pracma::repmat(TPT,N,1)  # LFTi=TPT (i=1..N)
  LST<- LFT-TD         # LSTi=LFTi-TDi (i=1..N)
  for (i in 1:N)       # Backward pass
  {
    TPTDSM=t(pracma::arrayfun("*",t(as.matrix(dsm)),LST)) #Calculate LST
    TPTDSM[is.nan(TPTDSM)*1]<- TPT         #Independent tasks' LFT = TPT
    TPTDSM[!TPTDSM]<-TPT              #Independent tasks' LFT = TPT
    LFT<-matrix(Rfast::rowMins(TPTDSM,value=TRUE))       #Calculate LFT
    LST[i]<-LFT[i]-TD[i]                       #LSTi=LFTi-TDi (i:=1..N)
  }
  output <- list(TPT=TPT,EST=EST,EFT=EFT,LST=LST,LFT=LFT)
  class(output)<-"TPT"
  return(output)
}
