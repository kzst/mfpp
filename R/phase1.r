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
#' @importFrom stats rbeta
phase1<- function(x,a=-0.1,b=0.30,pdftype="uniform"){
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop(
      "Package \"pracma\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package \"Matrix\" must be installed to use this function.",
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
        "phase1 works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  n=pracma::size(PDM,1)
  m=pracma::size(PDM,2)
  M=b-a
  PDMout<-PDM
  if ("uniform" %in% pdftype)
  {
    if (m>n){
      PDMout[,(n+1):m]=PDMout[,(n+1):m]+(M*pracma::rand(n,m-n)+a)*PDM[,(n+1):m]
      for (i in 1:n){                   #demands will be similar to the other demands
        for (j in ((n+1):m)){
          if ((PDM[i,j]>0) && (PDM[i,j]<=1) && (PDMout[i,j]>1))
            PDMout[i,j]=1                  # %Quality should not be greater than 1
        }
      }

    }
  }else{
    if ("beta" %in% pdftype){
      n <- pracma::size(PDM,1)
      m <- pracma::size(PDM,2)
      t <- (a+b)/2
      r1 <- a/t
      r2 <- b/t
      PDMout <- PDM
      if(b>a)
        alpha <- 6*(1-r1)/(r2-r1)
      beta <- 6*(r2-1)/(r2-r1)
      M <- b-a
      if (m>n){
        PDMout[,(n+1):m]=PDMout[,(n+1):m]+(M*matrix(rbeta(n*(m-n), alpha, beta), ncol=(m-n))+a)*PDM[,(n+1):m]
        for (i in 1:n){                   #demands will be similar to the other demands
          for (j in ((n+1):m)){
            if ((PDM[i,j]>0) && (PDM[i,j]<=1) && (PDMout[i,j]>1))
              PDMout[i,j]=1                  # %Quality should not be greater than 1
          }
        }

      }
    }else{
      warning("\n\nphase1 implemented only for 'uniform' and 'beta' distributions")
    }
  }
  class(PDMout)<-"PDM_matrix"
  if ("PDM_list" %in% class(x)){
    x$PDM<-PDMout
    output<-x
    class(output)<-"PDM_list"
    return(output)
  }else{
    return(PDMout)
  }
}




