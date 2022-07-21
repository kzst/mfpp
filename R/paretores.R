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
#' @importFrom utils tail

paretores<- function(DSM,TD,RD){
  output=list()
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
  if (!requireNamespace("nsga2R", quietly = TRUE)) {
    stop(
      "Package \"nsga2R\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (!requireNamespace("genalg", quietly = TRUE)) {
    stop(
      "Package \"genalg\" must be installed to use this function.",
      call. = FALSE
    )
  }
  DSM<-pracma::triu(round(DSM)) #DSM must be an upper triangular binary matrix
  T<-TD
  R<-RD
  N<-pracma::numel(TD)
  st<-tpt(DSM,TD)
  EST<-st$EST
  LST<-st$LST
  SST<-EST
  maxresfun<-function(SST){tpr(as.matrix(SST)[1:pracma::numel(TD)],DSM,
                               TD,as.matrix(RD))}
  if (dim(RD)[2]>1){ # Multi-objective case
    results<-nsga2R::nsga2R(fn=maxresfun,dim(RD)[1],dim(RD)[2],
                            lowerBounds = EST,upperBounds = LST)
    rd<-tail(results$objectives,n=1)
    colnames(rd)<-paste("R",1:ncol(rd),sep="_")
    rownames(rd)<-"TPR"
    SST<-t(as.matrix(tail(results$parameters,n=1)))

  }else{ # Single objective case
    results<-genalg::rbga(as.vector(EST),as.vector(LST),evalFunc=maxresfun)
    rd<-as.matrix(tail(results$best,n=1))
    colnames(rd)<-paste("R",1:ncol(rd),sep="_")
    rownames(rd)<-"TPR"
    SST<-t(as.matrix(tail(results$population,n=1)))
  }
  colnames(SST)<-"SST"
  output$RD<-rd
  output$SST<-SST
  return(output)
}

