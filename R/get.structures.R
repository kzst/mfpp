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
get.structures<- function(x,type=c("min","max","minimax","maximin","most")){
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
        "get.structures works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  N<-dim(PDM)[1]
  M<-dim(PDM)[2]
  if (N>M){
    stop(
      "number of rows must be less or equal than the columns",
      call. = FALSE
    )
  }else{
    output<-list()
    minPDM<-PDM
    minPDM[1:N,1:N]<-floor(minPDM[1:N,1:N])
    minPDM[(diag(minPDM)==0)*c(1:N),(diag(minPDM)==0)*c(1:N)]<-0
    class(minPDM)<-"PDM_matrix"
    maxPDM<-PDM
    maxPDM[1:N,1:N]<-ceiling(maxPDM[1:N,1:N])
    maxPDM[(diag(maxPDM)==0)*c(1:N),(diag(maxPDM)==0)*c(1:N)]<-0
    class(maxPDM)<-"PDM_matrix"
    mostPDM<-PDM
    mostPDM[1:N,1:N]<-round(mostPDM[1:N,1:N])
    mostPDM[(diag(mostPDM)==0)*c(1:N),(diag(mostPDM)==0)*c(1:N)]<-0
    class(mostPDM)<-"PDM_matrix"
    maximinPDM<-minPDM
    diag(maximinPDM)<-diag(maxPDM)
    maximinPDM[(diag(maximinPDM)==0)*c(1:N),(diag(maximinPDM)==0)*c(1:N)]<-0
    class(maximinPDM)<-"PDM_matrix"
    minimaxPDM<-maxPDM
    diag(minimaxPDM)<-diag(minPDM)
    minimaxPDM[(diag(minimaxPDM)==0)*c(1:N),(diag(minimaxPDM)==0)*c(1:N)]<-0
    class(minimaxPDM)<-"PDM_matrix"
    if ("min" %in% type){ # Calculate minimal structure
      minstruct<-list()
      minstruct$PDM<-minPDM
      if (methods::is(x,"PDM_list")){
        minstruct$w<-x$w
        minstruct$Rs<-x$Rs
        class(minstruct)<-"PDM_list"
        output$minstruct<-minstruct
      }else{
        output$minstruct<-minPDM
      }
    }
    if ("max" %in% type){ # Calculate maximal structure
      maxstruct<-list()
      maxstruct$PDM<-maxPDM
      if (methods::is(x,"PDM_list")){
        maxstruct$w<-x$w
        maxstruct$Rs<-x$Rs
        class(maxstruct)<-"PDM_list"
        output$maxstruct<-maxstruct
      }else{
        output$maxstruct<-maxPDM
      }
    }
    if ("most" %in% type){ # Calculate desired structure
      moststruct<-list()
      moststruct$PDM<-mostPDM
      if (methods::is(x,"PDM_list")){
        moststruct$w<-x$w
        moststruct$Rs<-x$Rs
        class(moststruct)<-"PDM_list"
        output$moststruct<-moststruct
      }else{
        output$moststruct<-mostPDM
      }
    }
    if ("minimax" %in% type){ # Calculate minimax structure
      minimaxstruct<-list()
      minimaxstruct$PDM<-minimaxPDM
      if (methods::is(x,"PDM_list")){
        minimaxstruct$w<-x$w
        minimaxstruct$Rs<-x$Rs
        class(minimaxstruct)<-"PDM_list"
        output$minimaxstruct<-minimaxstruct
      }else{
        output$minimaxstruct<-minimaxPDM
      }
    }
    if ("maximin" %in% type){ # Calculate maximin structure
      maximinstruct<-list()
      maximinstruct$PDM<-maximinPDM
      if (methods::is(x,"PDM_list")){
        maximinstruct$w<-x$w
        maximinstruct$Rs<-x$Rs
        class(maximinstruct)<-"PDM_list"
        output$maximinstruct<-maximinstruct
      }else{
        output$maximinstruct<-maximinPDM
      }
    }
    if (methods::is(x,"PDM_list")){
      class(output)<-"Set_PDM_list"
    }else{
      class(output)<-"Set_PDM_matrix"
    }
    return(output)
  }
}
