#-----------------------------------------------------------------------------#
#                                                                             #
#  MATRIX-BASED FLEXIBLE PROJECT PLANNING                                     #
#                bbn                                                             #
#  Written by: Zsolt T. Kosztyan, Aamir Saghir                                #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: May 2022                                                     #
#-----------------------------------------------------------------------------#
#' @export
truncpdm<- function(x){
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop(
      "Package \"pracma\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if ("PDM_list" %in% class(x)){
    PDM<-x$PDM
  }else{
    if (("PDM_matrix" %in% class(x))||("matrix" %in% class(x))||("array" %in% class(x))||("data.frame" %in% class(x))){
      PDM<-x
    }else{
      stop(
        "truncpdm works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  N<-dim(PDM)[1]
  M<-dim(PDM)[2]
  if (min(N,M)>0)
  {
    if (pracma::numel(which(diag(PDM)!=0,TRUE))>0){
      if (N>M){
        PDM<-PDM[(diag(PDM)!=0) * c(1:N),]
      }else{
        if (N<M){
          PDM<-PDM[(diag(PDM)!=0) * c(1:N),c((diag(PDM)!=0) * c(1:N),c((N+1):M))]
        }else{
          PDM<-PDM[(diag(PDM)!=0) * c(1:N),(diag(PDM)!=0) * c(1:N)]
        }
      }
    }else{
      PDM<-matrix(0,0,0)
      class(PDM)<-"PDM_matrix"
    }
  }else{
    PDM<-matrix(0,0,0)
    class(PDM)<-"PDM_matrix"
  }
  if ("PDM_list" %in% class(x)){
    x$PDM<-PDM
    output<-x
    class(output)<-"PDM_list"
    return(output)
  }else{
    class(PDM)<-"PDM_matrix"
    return(PDM)
  }
}
