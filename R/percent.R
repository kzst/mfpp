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
percent<- function(PDM,type=c("c","q","qd","r","s","t"),w=2,Rs=2,ratio=1){
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
  Const<-list()
  if ("PDM_list" %in% class(PDM)){
    Const$w<-PDM$w
    Const$Rs<-PDM$Rs
    w<-PDM$w
    Rs<-PDM$Rs
    PDM<-PDM$PDM
  }else{
    Const$w<-w
    Const$Rs<-Rs
  }
  Const$ratio<-ratio
  if ("c" %in% type)
  {
    DSMdiag <- matrix(ceiling(diag(PDM[,1:pracma::size(PDM,1)])))
    #All uncertain tasks/dependencies will be included

    dsmdiag <- matrix(floor(diag(PDM[,1:pracma::size(PDM,1)])))
    #All uncertain tasks/dependencies will be excluded

    C <- Rfast::rowMaxs(PDM[,
                    (pracma::size(PDM,1)+w+1):(pracma::size(PDM,1)+2*w)],
                    value=TRUE)
    c <- Rfast::rowMins(PDM[,
                    (pracma::size(PDM,1)+w+1):(pracma::size(PDM,1)+2*w)],
                    value=TRUE)
    TPCmax<- C%*%DSMdiag
    TPCmin <- c%*%dsmdiag
    if (TPCmax==TPCmin){
      Const$Cc <- as.numeric(TPCmin)
    }else{
      Const$Cc <- as.numeric(TPCmin+ratio*(TPCmax-TPCmin))
    }


  }
  if ("q" %in% type)
  {
    if (dim(PDM)[2]==dim(PDM)[1]+w*(3+Rs)) #There are QD
    {
      N <- pracma::size(PDM,1)             #Number of tasks
      PEM <- PDM[,1:N]                     #The original logic network
      DSM <- ceiling(PEM)                  #If all uncertainties are realized
      dsm <- floor(PEM)                    #If all uncertainties are ignored
      QD <- PDM[,(N+2*w+1):(N+3*w)]        #The quality domain
      Q <- matrix(Rfast::rowMaxs(QD, value=TRUE))  #The maximal quality level
      q <- matrix(Rfast::rowMins(QD, value=TRUE))  #The minimal quality level
      TPQmax <- tpq(DSM,PEM,Q)
      TPQmin <- tpq(dsm,PEM,q)
      Const$Cq <- TPQmin+ratio*(TPQmax-TPQmin)
    }
  }
  if ("qd" %in% type)
  {
    if (dim(PDM)[2]==dim(PDM)[1]+w*(3+Rs)) #There are QD
    {
      N <- pracma::size(PDM,1)             #Number of tasks
      PEM <- PDM[,1:N]                     #The original logic network
      DSM <- ceiling(PEM)                  #If all uncertainties are realized
      dsm <- floor(PEM)                    #If all uncertainties are ignored
      QD <- PDM[,(N+2*w+1):(N+3*w)]        #The quality domain
      Q <- matrix(Rfast::rowMaxs(QD, value=TRUE))  #The maximal quality level
      q <- matrix(Rfast::rowMins(QD, value=TRUE))  #The minimal quality level
      TPQmax <- tpq(DSM,PEM,QD,Q)
      TPQmin <- tpq(dsm,PEM,QD,q)
      Const$Cq <- TPQmin+ratio*(TPQmax-TPQmin)
    }
  }
  if ("r" %in% type)
  {
    if (dim(PDM)[2]==dim(PDM)[1]+w*(3+Rs)) #There are QD
    {
      DSM <- floor(pracma::triu(PDM[,1:pracma::size(PDM,1)],1))+diag(ceiling(diag(PDM)))  #If every
      #tasks will be included, however, every dependencies will be excluded
      dsm <- ceiling(pracma::triu(PDM[,1:pracma::size(PDM,1)],1))+diag(floor(diag(PDM)))    #If every
      #tasks will be excluded, however, every dependencies will be included
      rD <- PDM[,(pracma::size(PDM,1)+3*w+1):ncol(PDM)]
      R <- c()                           #Maximal values of resource demands
      r <- c()                           #Minimal values of resource demands
      if (w > 1)
        for (i in seq(1,pracma::size(rD,2),w)) {
          rmin <- matrix(Rfast::colMins(t(rD[,i:(i+w-1)]),value=TRUE))
          rmin <- stats::na.omit(rmin)
          rmax <- matrix(Rfast::colMaxs(t(rD[,i:(i+w-1)]),value=TRUE))
          rmax <- stats::na.omit(rmax)
          r <- cbind(r,rmin)
          R <- cbind(R,rmax)
        }  else {
          R <- rD
          r <- rD
        }
      T <- matrix(Rfast::rowMaxs(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))   #min R when max T
      t <- matrix(Rfast::rowMins(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))   #max R when min T
      EST <- tpt(DSM,t)[["EST"]]                                        #Optimization are within [EST,LST]
      LST <- tpt(DSM,t)[["LST"]]
      TPRmax=t(matrix(pmax(tpr(EST,DSM,t,as.matrix(R)),tpr(LST,DSM,t,as.matrix(R)))))
      if (ratio==1.0){
        CR=TPRmax
        colnames(CR)<-paste("R",1:ncol(CR),sep="_")
        rownames(CR)<-"TPR"
        Const$CR<-CR
        }  else {
          #calculation of TPRmin
          TPRmin=paretores(dsm,T,as.matrix(r))$RD
          Const$CR=TPRmin+ratio*(TPRmax-TPRmin)}

    }else{
      if (dim(PDM)[2]==dim(PDM)[1]+w*(2+Rs)) #There are no QD
      {
        DSM <- floor(pracma::triu(PDM[,1:pracma::size(PDM,1)],1))+diag(ceiling(diag(PDM)))  #If every
        #tasks will be included, however, every dependencies will be excluded
        dsm <- ceiling(pracma::triu(PDM[,1:pracma::size(PDM,1)],1))+diag(floor(diag(PDM)))    #If every
        #tasks will be excluded, however, every dependencies will be included
        rD <- PDM[,(pracma::size(PDM,1)+2*w+1):ncol(PDM)]
        R <- c()                           #Maximal values of resource demands
        r <- c()                           #Minimal values of resource demands
        if (w > 1)
          for (i in seq(1,pracma::size(rD,2),w)) {
            rmin <- matrix(Rfast::colMins(t(rD[,i:(i+w-1)]),value=TRUE))
            rmin <- stats::na.omit(rmin)
            rmax <- matrix(Rfast::colMaxs(t(rD[,i:(i+w-1)]),value=TRUE))
            rmax <- stats::na.omit(rmax)
            r <- cbind(r,rmin)
            R <- cbind(R,rmax)
          }
        else {
            R <- rD
            r <- rD
          }
        T <- matrix(Rfast::rowMaxs(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))   #min R when max T
        t <- matrix(Rfast::rowMins(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))   #max R when min T
        EST <- tpt(DSM,t)[["EST"]]                                        #Optimization are within [EST,LST]
        LST <- tpt(DSM,t)[["LST"]]
        TPRmax=t(matrix(pmax(tpr(EST,DSM,t,as.matrix(R)),tpr(LST,DSM,t,as.matrix(R)))))
        if (ratio==1.0){
          CR=TPRmax
          colnames(CR)<-paste("R",1:ncol(CR),sep="_")
          rownames(CR)<-"TPR"
          Const$CR<-CR
        }  else {
          #calculation of TPRmin
            TPRmin=paretores(dsm,T,as.matrix(r))$RD
            Const$CR<-TPRmin+ratio*(TPRmax-TPRmin)}
      }
    }
  }
  if ("s" %in% type)
  {
    PEM=PDM[,1:(pracma::size(PDM,1))]                         #N by N matrix of the logic domain
    TPSmax=maxscore_PEM(PEM,PEM,(pracma::ones(pracma::size(PEM,1)))-PEM)
    TPSmin=minscore_PEM(PEM,PEM,(pracma::ones(pracma::size(PEM,1)))-PEM)
    Const$Cs<-TPSmin+ratio*(TPSmax-TPSmin)

  }
  if ("t" %in% type)
  {
    DSM=ceiling(PDM[,1:pracma::size(PDM,1)])        #All uncertain tasks/dependencies will be included
    dsm=floor(PDM[,1:pracma::size(PDM,1)])          #All uncertain tasks/dependencies will be excluded
    T <- matrix(Rfast::rowMaxs(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))
    t <- matrix(Rfast::rowMins(PDM[,(pracma::size(PDM,1)+1):(pracma::size(PDM,1)+w)], value=TRUE))
    TPTmax=tpt(DSM,T)[[1]]
    TPTmin=tpt(dsm,t)[[1]]
    Const$Ct<-TPTmin+ratio*(TPTmax-TPTmin)

  }
  class(Const)<-"PDM_const"
  return(Const)
}







