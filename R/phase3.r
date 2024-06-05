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
phase3<- function(x,p=0.10,s=0.50,nW=0){
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
        "phase3 works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  n<-dim(PDM)[1]
  m<-dim(PDM)[2]
  PDMout=PDM
  PDMout[1:n,1:n]=pmax(pmin(pracma::ones(n),PDM[1:n,1:n]+s*pracma::triu(((pracma::rand(n)<p)*1)*pracma::rand(n))),pracma::zeros(n))

  if (m>n){                #occurances is generated with probability value p
    Z=pracma::zeros(n,(m-n))
    PDMout[,(n+1):m]=PDMout[,(n+1):m]+((PDM[,(n+1):m]==Z)*1)*pracma::rand(n,m-n)*pracma::repmat(colMeans(PDM[diag(PDM)!=0,(n+1):m]),n,1)  #is generated then
    for (i in 1:n){                   #demands will be similar to the other demands
      for (j in ((n+1):m)){
        if ((PDM[i,j]>0) && (PDM[i,j]<=1) && (PDMout[i,j]>1))
          PDMout[i,j]=1                  # %Quality should not be greater than 1
      }
    }
    PDMout[1:(n-nW),(n+1):m] <- PDM[1:(n-nW),(n+1):m] # Write back for regular tasks
    PDMout[diag(PDMout)==0,] <- 0        #Exluded task demands are also excluded
    PDMout[1:n, (diag(PDMout)==0)*c(1:n)]<- 0 #Exluded task demands are also excluded
  }
  class(PDMout)<-"PDM_matrix"
  if (methods::is(x,"PDM_list")){
    x$PDM<-PDMout
    output<-x
    class(output)<-"PDM_list"
    return(output)
  }else{
    return(PDMout)
  }
}




