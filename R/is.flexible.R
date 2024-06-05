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
is.flexible<- function(x){
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop(
      "Package \"pracma\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (methods::is(x,"PDM_list")){
    PDM<-x$PDM
  }else{
    if ((methods::is(x,"PDM_matrix"))||(methods::is(x,"matrix"))||(methods::is(x,"array"))||(methods::is(x,"data.frame"))){
      PDM<-x
    }else{
      stop(
        "is.flexible works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  N<-dim(PDM)[1]
  M<-dim(PDM)[2]
  if (min(N,M)>0)
  {
    if (N>M){
      return(FALSE) # Not a real PDM matrix
    }else{
      PEM<-PDM[1:N,1:N]
      if (pracma::numel(which(((PEM<1)&(PEM>0)),TRUE))>0){
        return(TRUE)
      }else{
        return(FALSE)
      }
    }
  }else{return(FALSE)}
}
