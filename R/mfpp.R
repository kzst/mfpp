#-----------------------------------------------------------------------------#
#                                                                             #
#  MATRIX-BASED FLEXIBLE PROJECT PLANNING                                     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan, Aamir Saghir                                #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: June 2024
# For reproducible results of papers
#----------------------------------------------------------------------------

########## Generation of project domain matrix (PDM) for flexible project planning###############

##### Function to generate a PDM ####
generatepdm<- function(N,ff,cf,mTD,mCD,mRD,w,nR,nW,scale=1.4,QD=FALSE,lst=FALSE)
{
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

  if(missing(scale))
  {scale <- 1.4}
  if(!QD) {
    cf <- cf+1
    PEM <- phase3(pracma::triu((pracma::triu(pmin(pracma::ones(N)/
                                                    pmax(pracma::repmat((1-cf):(N-cf),N,1)^(scale)-
                                                           (pracma::repmat(matrix(0:(N-1)),1,N)^(scale)),
                                                         pracma::ones(N)),pracma::ones(N)),1) > pracma::rand(N))*1,1)+
                    pracma::eye(N),ff,-0.5)   # Generate PEM matrix
    nTD <- w                          #Width of TD = number of modes
    nCD <- w                          #Width of CD = number of modes
    nRD <- w*nR                       #Width of RD = number of modes x number of
    #resources
    TD <- pracma::rand(N,nTD)*mTD             #Generate time domain
    CD <- pracma::rand(N,nCD)*mCD             #Generate cost domain
    rD <- pracma::rand(N,nRD)*mRD             #Generate resource domain
    pem <- pracma::zeros(N+nW)
    pem[1:N,1:N]<- PEM
    td <- pracma::zeros(N+nW,nTD)
    cd <- pracma::zeros(N+nW,nCD)
    rd <- pracma::zeros(N+nW,nRD)
    if (w==2)                       #In case of CTCTP the columns will be sorted
    {
      TD <- cbind(matrix(Rfast::rowMaxs(TD, value=TRUE))-
                    matrix(Rfast::rowMaxs(TD, value=TRUE))*
                    pracma::rand(N,1)*0.2,matrix(Rfast::rowMaxs(TD, value=TRUE)))
      CD <- cbind(matrix(Rfast::rowMaxs(CD, value=TRUE))-
                    matrix(Rfast::rowMaxs(CD, value=TRUE))*
                    pracma::rand(N,1)*0.2,matrix(Rfast::rowMaxs(TD, value=TRUE)))
      RD <- c()
      for (i in seq(1,nRD,2)){
        rmax <- matrix(Rfast::colMaxs(t(rD[,i:(i+1)]),value=TRUE))
        rmax <- na.omit(rmax)
        rmin <- rmax-rmax*pracma::rand(N,1)*0.2
        rmin <- na.omit(rmin)
        RD <- cbind(RD,rmin,rmax) }
    } else {
      RD <- rD
    }
    td[1:N,1:nTD] <- TD
    cd[1:N,1:nCD] <- CD
    rd[1:N,1:nRD] <- RD
    PDM <- cbind(pem,td,cd,rd)

    Rs<-nR
    rownames(PDM)<-paste("a",1:nrow(PDM),sep="_")
    if (Rs>0){
      colnames(PDM)<-c(paste("a",1:nrow(PDM),sep='_'),
                       paste("t",1:w,sep='_'),paste("c",1:w,sep='_'),
                       paste(paste("r",1:w,sep='_'),rep(1:Rs,each=w),sep='.'))
    }else{
      colnames(PDM)<-c(paste("a",1:nrow(PDM),sep='_'),
                       paste("t",1:w,sep='_'),paste("c",1:w,sep='_'))
    }

    class(PDM)<-"PDM_matrix"
    output$PDM<- PDM
    output$w <- w
  } else {
    cf=cf+1
    PEM=phase3(pracma::triu((pracma::triu(pmin(pracma::ones(N)/
                                                 pmax(pracma::repmat((1-cf):(N-cf),N,1)^(scale)-
                                                        (pracma::repmat(matrix(0:(N-1)),1,N)^(scale)),pracma::ones(N)),
                                               pracma::ones(N)),1) > pracma::rand(N))*1,1)+pracma::eye(N),ff,-0.5)
    # Generate PEM matrix
    nTD=w                          #Width of TD = number of modes
    nCD=w                          #Width of CD = number of modes
    nQD=w                          #Width of QD = number of modes
    nRD=w*nR                       #Width of RD = number of modes
    #x number of resources
    TD=pracma::rand(N,nTD)*mTD             #Generate time domain
    CD=pracma::rand(N,nCD)*mCD             #Generate cost domain
    QD=pracma::rand(N,nQD)                 #Generate quality domain
    rD=pracma::rand(N,nRD)*mRD             #Generate resource domain
    pem=pracma::zeros(N+nW)
    pem[1:N,1:N]=PEM
    td=pracma::zeros(N+nW,nTD)
    cd=pracma::zeros(N+nW,nCD)
    qd=pracma::zeros(N+nW,nQD)
    rd=pracma::zeros(N+nW,nRD)
    if (w==2)                       #In case of CTCTP the columns will be sorted
    {
      TD=cbind(matrix(Rfast::rowMaxs(TD, value=TRUE))-
                 matrix(Rfast::rowMaxs(TD, value=TRUE))*
                 pracma::rand(N,1)*0.2,matrix(Rfast::rowMaxs(TD, value=TRUE)))
      CD=cbind(matrix(Rfast::rowMaxs(CD, value=TRUE))-
                 matrix(Rfast::rowMaxs(CD, value=TRUE))*
                 pracma::rand(N,1)*0.2,matrix(Rfast::rowMaxs(CD, value=TRUE)))
      QD=cbind(matrix(Rfast::rowMaxs(QD, value=TRUE))-
                 matrix(Rfast::rowMaxs(QD, value=TRUE))*
                 pracma::rand(N,1)*0.2,matrix(Rfast::rowMaxs(QD, value=TRUE)))
      RD=c()
      for (i in seq(1,nRD,2)){
        rmax=matrix(Rfast::colMaxs(t(rD[,i:(i+1)]),value=TRUE))
        rmax=na.omit(rmax)
        rmin=rmax-rmax*pracma::rand(N,1)*0.2
        rmin=na.omit(rmin)
        RD=cbind(RD,rmin,rmax) }
    } else {
      RD=rD
    }
    td[1:N,1:nTD]=TD
    cd[1:N,1:nCD]=CD
    qd[1:N,1:nQD]=QD
    rd[1:N,1:nRD]=RD
    PDM=cbind(pem,td,cd,qd,rd)
    Rs<-nR
    rownames(PDM)<-paste("a",1:nrow(PDM),sep="_")
    if (Rs>0){
      colnames(PDM)<-c(paste("a",1:nrow(PDM),sep='_'),
                       paste("t",1:w,sep='_'),paste("c",1:w,sep='_'),
                       paste("q",1:w,sep='_'),
                       paste(paste("r",1:w,sep='_'),rep(1:Rs,each=w),sep='.'))
    }else{
      colnames(PDM)<-c(paste("a",1:nrow(PDM),sep='_'),
                       paste("t",1:w,sep='_'),paste("c",1:w,sep='_'),
                       paste("q",1:w,sep='_'))
    }
    class(PDM)<-"PDM_matrix"
    output$PDM<-PDM
    output$w<-w
  }
  class(PDM)<-"PDM_matrix"
  if (lst==FALSE){
    return(PDM)
  }else{
    output<-list()
    output$PDM<-PDM
    output$w<-w
    output$Rs<-nR
    class(output)<-"PDM_list"
    return(output)
  }
}

#### Function to calculate minimal/maximal/most likely project structures ####
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

#### Function to check the  flexibility of PDM matrix####

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

####Function to drop excluded tasks in PDM####
truncpdm<- function(x){
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
  if (methods::is(x,"PDM_list")){
    x$PDM<-PDM
    output<-x
    class(output)<-"PDM_list"
    return(output)
  }else{
    class(PDM)<-"PDM_matrix"
    return(PDM)
  }
}

############Calculation of project scores######################

####Function to calculate maximal score value (PMAX) of possible project scenarios####

maxscore_PEM<- function(PEM,P=PEM, Q=1-PEM)
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
  score=1
  p=diag(P)
  q=diag(Q)
  pem=diag(PEM)
  N=pracma::numel(pem)
  pqmax=Rfast::rowMaxs(matrix(c(p,q),ncol=2),value=TRUE)
  if (N>0)
    #The score of the project scenario is the geometric mean of maximum
    score=prod(matrix(c(p[pem==1], q[pem==0], pqmax[pem>0 & pem<1])))^{1/N}
  return(score)
}

####Function to calculate minimal score value of possible project scenarios####

minscore_PEM<- function(PEM,P=PEM, Q=1-PEM)
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

  score=1
  N=0
  p=diag(P)
  q=diag(Q)
  pem=diag(PEM)
  N=pracma::numel(pem)
  pqmin=Rfast::rowMins(matrix(c(p,q),ncol=2),value=TRUE)
  # for(i in 1:N) {
  #if ((p[i]>0) & (p[i]<1))
  # N=N+1
  #if (pem[i]==1)
  # score=score*p[i]
  #if (pem[i]==0)
  # score=score*q[i]
  #if (pem[i]<1 & pem[i]>0)
  # score=score*min(p[i],q[i])}

  if (N>0)           #The score of the project scenario is the geometric mean of maximum
    score=prod(matrix(c(p[pem==1], q[pem==0], pqmin[pem>0 & pem<1])))^{1/N}

  return(score)
}

#### Function to calculate Pareto-optimal resource allocation####
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

####Function to calculate desired project completion characteristic of a project structure ####

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
  if (methods::is(PDM,"PDM_list")){
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
          rmin <- na.omit(rmin)
          rmax <- matrix(Rfast::colMaxs(t(rD[,i:(i+w-1)]),value=TRUE))
          rmax <- na.omit(rmax)
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
            rmin <- na.omit(rmin)
            rmax <- matrix(Rfast::colMaxs(t(rD[,i:(i+w-1)]),value=TRUE))
            rmax <- na.omit(rmax)
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

######### Sensitivity analysis of PDM #######################

#### Function to simulate estimation uncertainty####
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
  if (methods::is(x,"PDM_list")){
    PDM<-x$PDM
  }else{
    if ((methods::is(x,"PDM_matrix"))||(methods::is(x,"matrix"))||(methods::is(x,"array"))||(methods::is(x,"data.frame"))){
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
  if (methods::is(x,"PDM_list")){
    x$PDM<-PDMout
    output<-x
    class(output)<-"PDM_list"
    return(output)
  }else{
    return(PDMout)
  }
}

####Function to simulate shock effects#####
phase2<- function(x,p=0.1,s=5.0){
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
        "truncpdm works only on matix, PDM_matrix, and PDM_list.",
        call. = FALSE
      )
    }
  }
  class(PDM)<-"PDM_matrix"
  n=pracma::size(PDM,1)
  m=pracma::size(PDM,2)
  PDMout=PDM
  if (m>n){
    PDMout[,(n+1):m]=PDMout[,(n+1):m]+(pracma::rand(n,m-n)<p)*(PDM[,(n+1):m]*s)
    for (i in 1:n){                   #demands will be similar to the other demands
      for (j in ((n+1):m)){
        if ((PDM[i,j]>0) && (PDM[i,j]<=1) && (PDMout[i,j]>1))
          PDMout[i,j]=1                  # %Quality should not be greater than 1
      }
    }

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

####Function to simulate the effects of the change of customer claims####
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


##### Plot function for Matrix-Based Flexible Project Planning####
#' @export
#' @importFrom graphics par legend barplot
plot.PDM_matrix <- function(x,w=NULL,Rs=NULL,
                            type=c("orig","max","min","maximin","minimax","most","const"),
                            main=NULL,col=NULL,
                            ...){
  if (methods::is(x,"PDM_matrix")){
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop(
        "Package \"igraph\" must be installed to use this function.",
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
    oldpar<-par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(1,1))
    PDM<-x
    class(PDM)<-"PDM_matrix"
    N<-pracma::size(PDM,1)
    if (is.null(rownames(PDM))) rownames(PDM)<-paste("a",1:N,sep="_")
    M<-pracma::size(PDM,2)
    if (N>M){
      stop(
        "number of rows must be less or equal than the columns",
        call. = FALSE
      )
    }else{
      if (is.flexible(PDM)){
        pdm<-truncpdm(PDM)
        n<-pracma::size(pdm,1)
        c<-which((diag(pdm)<1)&(diag(pdm)>0),TRUE)
        diag(pdm)<-0
        g<-igraph::graph.adjacency(pdm[1:n,1:n],weighted = TRUE)
        igraph::E(g)$color="black"
        if (pracma::numel(which(igraph::E(g)$weight<1,TRUE))){
          igraph::E(g)[igraph::E(g)$weight<1]$color="grey"
        }
        igraph::V(g)$color="green"
        if (pracma::numel(c)>0){
          igraph::V(g)[c]$color="yellow"
        }
        if ("orig" %in% type){
          if (!is.null(main)){
            plot(g,main=main,
                 layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,...)
          }else{
            plot(g,main="Original Logic Network",
                 layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,...)
          }
          legend(
            "topleft",
            legend = c("mandatory", "supplementary"),
            pt.bg  = c("green", "yellow"),
            pch    = 22,
            cex    = 1,
            bty    = "n",
            title  = "Tasks"
          )
          legend(
            "bottomleft",
            legend = c("fixed", "flexible"),
            col  = c("black", "grey"),
            pch    = 45,
            cex    = 1,
            bty    = "n",
            title  = "Dependencies"
          )
        }
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

        minpdm<-truncpdm(minPDM)
        n<-pracma::size(minpdm,1)
        m<-pracma::size(minpdm,2)
        c<-NULL
        if (!is.null(w)) { # Number of completion mode is specified
          if (m>=(n+w)){ # There are a Task Domain
            TPT<-tpt(minpdm[1:n,1:n],Rfast::rowMins(minpdm[,(n+1):(n+w)]))
            c<-which(as.vector(TPT$EFT)==as.vector(TPT$LFT),TRUE)
          }
        }
        diag(minpdm)<-0
        g<-igraph::graph.adjacency(minpdm[1:n,1:n],weighted=TRUE)
        if (!is.null(w)) igraph::V(g)$weight<-Rfast::rowMins(minpdm[,(n+1):(n+w)])
        igraph::V(g)$color="green"
        if (!is.null(c)) igraph::V(g)[c]$color="red"
        if ("min" %in% type){
          if (!is.null(c)){
            plot(g,main="Minimal Structure",layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,
                 vertex.label=paste("d",igraph::V(g)$weight,sep="="),...)

            legend(
              "topleft",
              legend = c("critical", "non-critical"),
              pt.bg  = c("red", "green"),
              pch    = 22,
              cex    = 1,
              bty    = "n",
              title  = "Tasks"
            )
          }else{
            plot(g,main="Minimal Structure",layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,...)

          }
        }
        maxpdm<-truncpdm(maxPDM)
        n<-pracma::size(maxpdm,1)
        m<-pracma::size(maxpdm,2)
        c<-NULL
        if (!is.null(w)){ # Number of completion mode is specified
          if (m>=(n+w)){ # There are a Task Domain
            TPT<-tpt(maxpdm[1:n,1:n],Rfast::rowMaxs(maxpdm[,(n+1):(n+w)]))
            c<-which(as.vector(TPT$EFT)==as.vector(TPT$LFT),TRUE)
          }
        }
        diag(maxpdm)<-0
        g<-igraph::graph.adjacency(maxpdm[1:n,1:n],weighted=TRUE)
        if (!is.null(w)) igraph::V(g)$weight<-Rfast::rowMaxs(maxpdm[,(n+1):(n+w)])
        igraph::V(g)$color="green"
        if (!is.null(c)) igraph::V(g)[c]$color="red"

        if ("max" %in% type){
          if (!is.null(c)){
            plot(g,main="Maximal Structure",layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,
                 vertex.label=paste("d",
                                    igraph::V(g)$weight,sep="="),...)

            legend(
              "topleft",
              legend = c("critical", "non-critical"),
              pt.bg  = c("red", "green"),
              pch    = 22,
              cex    = 1,
              bty    = "n",
              title  = "Tasks"
            )
          }else{
            plot(g,main="Maximal Structure",layout=igraph::layout_as_tree,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,...)

          }
        }
        minimaxpdm<-truncpdm(minimaxPDM)
        diag(minimaxpdm)<-0
        n<-pracma::size(minimaxpdm,1)
        if ("minimax" %in% type){
          plot(igraph::graph.adjacency(minimaxpdm[1:n,1:n]),main="Minimax Structure",
               layout=igraph::layout_as_tree,vertex.shape="crectangle",vertexlabel.dist=2.5,
               vertex.color="green",...)
        }
        maximinpdm<-truncpdm(maximinPDM)
        diag(maximinpdm)<-0
        n<-pracma::size(maximinpdm,1)
        if ("maximin" %in% type){
          plot(igraph::graph.adjacency(maximinpdm[1:n,1:n]),
               main="Maximin Structure",
               layout=igraph::layout_as_tree,
               vertex.shape="crectangle",vertex.color="green",vertexlabel.dist=2.5,...)
        }
        mostpdm<-truncpdm(mostPDM)
        diag(mostpdm)<-0
        n<-pracma::size(mostpdm,1)
        if ("most" %in% type){
          plot(igraph::graph.adjacency(mostpdm[1:n,1:n]),
               main="Most-likely/Most-desired Structure",
               layout=igraph::layout_as_tree,
               vertex.shape="crectangle",vertex.color="green",
               vertexlabel.dist=2.5,...)
        }
      }else{ # For non flexible structures
        pdm<-truncpdm(PDM)
        c<-which((diag(pdm)<1)&(diag(pdm)>0),TRUE)
        diag(pdm)<-0
        n<-pracma::size(pdm,1)
        g<-igraph::graph.adjacency(pdm[1:n,1:n],weighted = TRUE)
        igraph::E(g)$color="black"
        if (pracma::numel(which(igraph::E(g)$weight<1,TRUE))){
          igraph::E(g)[igraph::E(g)$weight<1]$color="grey"
        }
        igraph::V(g)$color="green"
        if (pracma::numel(c)>0){
          igraph::V(g)[c]$color="yellow"
        }
        if (is.null(rownames(PDM))) rownames(PDM)<-paste("a",1:nrow(PDM),sep="_")
        igraph::V(g)$names<-rownames(PDM)
        if ("orig" %in% type){
          if (!is.null(main)){
            plot(g,main=main,
                 layout=igraph::layout_as_tree,vertex.shape="crectangle",
                 vertex.label=igraph::V(g)$names,vertexlabel.dist=2.5,...)
          }else{
            plot(g,main="Logic Network",
                 layout=igraph::layout_as_tree,vertex.label=igraph::V(g)$names,
                 vertex.shape="crectangle",vertexlabel.dist=2.5,...)
          }



        }
      }
    }
    if ("const" %in% type){
      type<-c("c","q","r","s","t")
      if (is.null(w)||is.null(Rs)){
        type<-"s"
      }else{
        if (Rs==0){
          type<-c("c","q","s","t")
        }
      }
      minCONST<-percent(PDM,type=type,w=w,Rs=Rs,ratio=0)
      maxCONST<-percent(PDM,type=type,w=w,Rs=Rs,ratio=1.0)
      n<-length(minCONST)-3
      if (n>0){
        oldpar<-par(no.readonly = TRUE)
        on.exit(par(oldpar))
        par(mfrow=c(1,n))
      }
      if (!is.null(minCONST$Ct)&&!is.null(maxCONST$Ct))
        barplot(cbind(minCONST$Ct,maxCONST$Ct),
                names.arg = c("min_C_t","max_C_t"),
                ylab = "TPT",main = "Duration constraints",
                col=col)
      if (!is.null(minCONST$Cc)&&!is.null(maxCONST$Cc))
        barplot(cbind(minCONST$Cc,maxCONST$Cc),
                names.arg = c("min_C_c","max_C_c"),
                ylab = "TPC",main = "Cost constraints",
                col=col)
      if (!is.null(minCONST$Cq)&&!is.null(maxCONST$Cq))
        barplot(cbind(minCONST$Cq,maxCONST$Cq),
                names.arg = c("min_C_q","max_C_q"),
                ylab = "TPQ",main = "Quality constraints",
                col=col)
      if (!is.null(minCONST$Cs)&&!is.null(maxCONST$Cs))
        barplot(cbind(minCONST$Cs,maxCONST$Cs),
                names.arg = c("min_C_s","max_C_s"),
                ylab = "TPS",main = "Scope/score constraints",
                col=col)
      if (!is.null(minCONST$CR)&&!is.null(maxCONST$CR))
        barplot(cbind(minCONST$CR,maxCONST$CR),
                names.arg = c(paste("min",colnames(minCONST$CR),sep="_"),
                              paste("max",colnames(maxCONST$CR),sep="_")),
                ylab = "TPR",main = "Resource constraints",
                col=col)
    }
  }else{
    plot(x,...)
  }
}

#' @export
plot.PDM_list <- function(x,
                          type=c("orig","max","min","maximin","minimax","most","const"),
                          main=NULL,col=NULL,
                          ...){
  if (methods::is(x,"PDM_list")){
    plot.PDM_matrix(x=x$PDM,w=x$w,Rs=x$Rs,
                    type=type,main=main,col=col,
                    ...)
  }else{
    plot(x,...)
  }
}

#' @export
plot.Set_PDM_matrix <- function(x,w=NULL,Rs=NULL,
                                type=c("orig","max","min",
                                       "maximin","minimax","most","const"),
                                col=NULL,
                                ...){
  if (methods::is(x,"Set_PDM_matrix")){
    if (!is.null(x$minstruct))
      plot.PDM_matrix(x=x$minstruct,w=w,Rs=Rs,
                      type=type,main="Minimal Structure",col=col,
                      ...)
    if (!is.null(x$maxstruct))
      plot.PDM_matrix(x=x$maxstruct,w=w,Rs=Rs,
                      type=type,main="Maximal Structure",col=col,
                      ...)
    if (!is.null(x$minimaxstruct))
      plot.PDM_matrix(x=x$minimaxstruct,w=w,Rs=Rs,
                      type=type,main="Minimax Structure",col=col,
                      ...)
    if (!is.null(x$maximinstruct))
      plot.PDM_matrix(x=x$maximinstruct,w=w,Rs=Rs,
                      type=type,main="Maximin Structure",col=col,
                      ...)
    if (!is.null(x$moststruct))
      plot.PDM_matrix(x=x$moststruct,w=w,Rs=Rs,
                      type=type,main="Most-likely/Most-desired Structure",
                      col=col,
                      ...)
  }else{
    plot(x,...)
  }
}

#' @export
plot.Set_PDM_list <- function(x,type=c("orig","max",
                                       "min","maximin",
                                       "minimax","most","const"),
                              col=NULL,
                              ...){
  if (methods::is(x,"Set_PDM_list")){
    if (!is.null(x$minstruct))
      plot.PDM_list(x=x$minstruct,
                    type=type,main="Minimal Structure",col=col,
                    ...)
    if (!is.null(x$maxstruct))
      plot.PDM_list(x=x$maxstruct,
                    type=type,main="Maximal Structure",col=col,
                    ...)
    if (!is.null(x$minimaxstruct))
      plot.PDM_list(x=x$minimaxstruct,main="Minimax Structure",
                    col=col,
                    type=type,
                    ...)
    if (!is.null(x$maximinstruct))
      plot.PDM_list(x=x$maximinstruct,
                    type=type,main="Maximin Structure",col=col,
                    ...)
    if (!is.null(x$moststruct))
      plot.PDM_list(x=x$moststruct,
                    type=type,main="Most-likely/Most-desired Structure",
                    col=col,
                    ...)
  }else{
    plot(x,...)
  }
}

#' @export
plot.TPT <- function(x,sched="E",...){
  if (methods::is(x,"TPT")){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "Package \"ggplot2\" must be installed to use this function.",
        call. = FALSE
      )
    }
    if (!requireNamespace("reshape2", quietly = TRUE)) {
      stop(
        "Package \"reshape2\" must be installed to use this function.",
        call. = FALSE
      )
    }
    ST<-as.matrix(x$EST)
    FT<-as.matrix(x$EFT)
    Crit<-as.matrix(x$EST==x$LST)
    if ("L" %in% sched) {
      ST<-as.matrix(x$LST)
      FT<-as.matrix(x$LFT)
    }
    if ("S" %in% sched) {
      ST<-as.matrix(x$SST)
      FT<-as.matrix(x$SFT)
    }
    if (is.null(rownames(ST))) rownames(ST)<-paste("a",1:nrow(ST),sep="_")
    value<-name<-is.critical<-start.date<-end.date<-NULL
    df<-data.frame(name=factor(rownames(ST),levels=rownames(ST)),start.date=ST,end.date=FT,is.critical=Crit)
    colnames(df)<-c("name","start.date","end.date","is.critical")
    mdf<-reshape2::melt(df, measure.vars = c("start.date", "end.date"))
    ggplot2::ggplot(mdf, ggplot2::aes(value, name, colour = is.critical)) +
      ggplot2::geom_line(size = 6) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL)
  }else{
    plot(x,...)
  }
}

####Summary function to print PDM constraints, matrices, lists, sets, collections####
summary.PDM_const <- function(object, digits =  getOption("digits"), ...) {
  if (methods::is(object,"PDM_const")){
    cat("\nSummary of the PDM constraints structure:\n")
    if (!is.null(object$Ct)) cat("\nTime constraint (Ct): ",
                                 round(object$Ct,digits = digits))
    if (!is.null(object$Cc)) cat("\nCost constraint (Cc): ",
                                 round(object$Cc,digits =digits))
    if (!is.null(object$Cs)) cat("\nScore/scope constraint (Cs): ",
                                 round(object$Cs,digits =digits))
    if (!is.null(object$Cq)) cat("\nQuality constraint (Cq): ",
                                 round(object$Cq,digits =digits))
    if (!is.null(object$CR)) {
      cat("\nResource constraint(s) (CR):\n")
      round(object$CR,digits=digits)
    }
  }else{
    summary(object,digits=digits,...)
  }
}

#' @export
summary.PDM_matrix <- function(object, digits =  getOption("digits"),
                               w=getOption("w"),
                               Rs=getOption("Rs"), ...) {
  if (methods::is(object,"PDM_matrix")){
    cat("\nsummary PDM matrix:\n")
    print(object,digits=digits)
    if ((!is.null(w))||(!is.null(Rs))){
      if (Rs>0){
        maxCONST<-percent(object,type=c("c","q","r","s","t"),w=w,Rs=Rs,
                          ratio=1)
      }else{
        maxCONST<-percent(object,type=c("c","q","s","t"),w=w,Rs=Rs,
                          ratio=1)
      }
      minCONST<-percent(object,type=c("c","q","s","t"),w=w,Rs=Rs,
                        ratio=0)
      cat("\nMinimal constraints:\n")
      summary.PDM_const(minCONST,digits=digits)
      cat("\n\nMaximal constraints:\n")
      summary.PDM_const(maxCONST,digits=digits)
    }else{
      maxCONST<-percent(object,type=c("s"),ratio=1)
      minCONST<-percent(object,type=c("s"),ratio=0)
      cat("\nMinimal constraints:\n")
      summary.PDM_const(minCONST,digits=digits)
      cat("\n\nMaximal constraints:\n")
      summary.PDM_const(maxCONST,digits=digits)
      cat("\n\n")
      warning("Number of completion modes (w), and ",
              "number of resources (Rs) should be specified ",
              "to calculate constraints of demands.")
    }
  }else{
    summary(object,digits=digits,w=w,Rs=Rs,...)
  }
}

#' @export
summary.PDM_list <- function(object, digits =  getOption("digits"), ...) {
  if (methods::is(object,"PDM_list")){
    cat("\nsummary PDM list:\n")
    if (!is.null(object$w)) cat("\nNumber of completion modes (w): ",
                                object$w)
    if (!is.null(object$Rs)) cat("\nNumber of resources (Rs): ",
                                 object$Rs)
    summary.PDM_matrix(object$PDM,digits=digits,w=object$w,Rs=object$Rs)
  }else{
    summary(object,digits=digits,...)
  }
}

#' @export
summary.Set_PDM_matrix <- function(object, digits =  getOption("digits"),
                                   w=getOption("w"),
                                   Rs=getOption("Rs"), ...) {
  if (methods::is(object,"Set_PDM_matrix")){
    cat("\nSummary of main structures:\n")
    if (!is.null(object$minstruct)) {
      cat("\n\n\nSummary of minimal structure:\n")
      summary.PDM_matrix(object$minstruct,digits = digits,w=w,Rs=Rs)
    }
    if (!is.null(object$maxstruct)) {
      cat("\n\n\nSummary of maximal structure:\n")
      summary.PDM_matrix(object$maxstruct,digits = digits,w=w,Rs=Rs)
    }
    if (!is.null(object$minimaxstruct)) {
      cat("\n\n\nSummary of minimax structure:\n")
      summary.PDM_matrix(object$minimaxstruct,digits = digits,w=w,Rs=Rs)
    }
    if (!is.null(object$maximinstruct)) {
      cat("\n\n\nSummary of maximin structure:\n")
      summary.PDM_matrix(object$maximinstruct,digits = digits,w=w,Rs=Rs)
    }
    if (!is.null(object$moststruct)) {
      cat("\n\n\nSummary of most likely/most desired structure:\n")
      summary.PDM_matrix(object$moststruct,digits = digits,w=w,Rs=Rs)
    }
  }else{
    summary(object,digits=digits,w=w,Rs=Rs,...)
  }
}

#' @export
summary.Set_PDM_list <- function(object, digits =  getOption("digits"), ...) {
  if (methods::is(object,"Set_PDM_list")){
    cat("\nSummary of main structures:\n")
    if (!is.null(object$minstruct)) {
      cat("\n\n\nSummary of minimal structure:\n")
      summary.PDM_list(object$minstruct,digits = digits)
    }
    if (!is.null(object$maxstruct)) {
      cat("\n\n\nSummary of maximal structure:\n")
      summary.PDM_list(object$maxstruct,digits = digits)
    }
    if (!is.null(object$minimaxstruct)) {
      cat("\n\n\nSummary of minimax structure:\n")
      summary.PDM_list(object$minimaxstruct,digits = digits)
    }
    if (!is.null(object$maximinstruct)) {
      cat("\n\n\nSummary of maximin structure:\n")
      summary.PDM_list(object$maximinstruct,digits = digits)
    }
    if (!is.null(object$moststruct)) {
      cat("\n\n\nSummary of most likely/most desired structure:\n")
      summary.PDM_list(object$moststruct,digits = digits)
    }
  }else{
    summary(object,digits=digits,...)
  }
}

#' @export
summary.Collection_PDM <- function(object, digits =  getOption("digits"), ...) {
  if (methods::is(object,"Collection_PDM")){
    cat("\n\n\nSummary of PDM collection:\n")
    cat("\nNumber of projects: ",length(object))
    cat("\nList of projects: ")
    df<-data.frame(Project_name=names(object))
    df$w<-1
    df$Rs<-0
    for (i in (1:length(object))){
      df[i,"w"]<-as.numeric(object[[i]]$PDM_list$w)
      df[i,"Rs"]<-as.numeric(object[[i]]$PDM_list$Rs)

      cat("\n Project name: ",names(object)[i])
      cat(", w: ",df[i,"w"])
      cat(", Rs: ",df[i,"Rs"])
    }
    # invisible(df)
    return(invisible(df))
  }else{
    summary(object,digits=digits,...)
  }
  invisible()
}

#' @export
summary.TPT <- function(object, digits =  getOption("digits"),...){
  if (methods::is(object,"TPT")){
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop(
        "Package \"knitr\" must be installed to use this function.",
        call. = FALSE
      )
    }
    cat("\n\n Table of schedule\n")
    TPT<-object
    if (is.null(rownames(TPT$EST))) rownames(TPT$EST)<-paste("a",1:nrow(TPT$EST),sep="_")
    df<-data.frame(duration=TPT$EFT-TPT$EST,
                   EST=TPT$EST,EFT=TPT$EFT,LST=TPT$LST,LFT=TPT$LFT,
                   TF=TPT$LST-TPT$EST,SST=TPT$SST,SFT=TPT$SFT,
                   SF=TPT$LST-TPT$SST,Crit=TPT$EST==TPT$LST)
    colnames(df)<-c("Dur","EST","EFT","LST","LFT","TF","SST","SFT","SF","Is.Crit")
    print.data.frame(df,digits=digits)
    return(invisible(df))
  }else{
    summary(object,digits=digits,...)
  }
  invisible(df)
}

######## Calculation of Project demands##################

####Function of Cost demands of a project#####
tpc <- function(DSM,CD){
  TPC<-0
  for (i in 1:length(CD)){
    TPC=TPC+as.numeric(DSM[i,i])*as.numeric(CD[i])
  }
  return(as.numeric(TPC))
}

####Function to calcualte Total Project Quality for a project structure#####

tpq<- function(DSM,PEM, q, QD=NULL)
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

  TPQ <- 0  # Total Project Quality
  TPS <- 0  # Total Project Score (additive scores are assumed)

  for (i in 1:ncol(DSM))
  {
    if (DSM[i,i]>0)
    {TPS <- TPS+PEM[i,i]}
    if (TPQ==0)
    {TPQ=1}

    TPQ <- TPQ*(q[i]^PEM[i,i])
  }
  if (TPS>0)
    TPQ <- maxscore_PEM(DSM,PEM,(pracma::ones(pracma::size(PEM,2))-PEM))*(TPQ^{1/TPS})
  if (is.null(QD)){
    output <- TPQ
  }else{
    ## CONT
    pem <- matrix(diag(PEM))
    dsm <- matrix(diag(DSM))
    TPQ <- 0
    if (sum(Rfast::rowMaxs(QD[pem>0,], value = TRUE))>0)
      TPQ <- sum(q[dsm>0])/sum(Rfast::rowMaxs(QD[pem>0,], value = TRUE))
    output <-TPQ
  }
  return(output)
}

####Function to calculate maximum resource demands of a project####
tpr<- function(SST,DSM,TD,RD,res.graph=FALSE)
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
  SST[is.na(SST)]<-0
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
            {
              SST[j]<- SFT[i]
              SFT[j]<- SST[j]+TD[j]
            }
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
  if (res.graph==TRUE){
    TPT<-tpt(DSM,TD,SST)$TPT
    Width=as.vector(matrix(0,1,length(BP)))
    if (length(BP)>1)
      for (i in 2:length(BP)){
        Width[i-1]<-BP[i]-BP[i-1]
      }
    Width[length(BP)]<-TPT-BP[length(BP)]
    m<-pracma::size(RESFUNC,2)
    oldpar<-par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow=c(m,1))
    for (i in c(1:m)){
      barplot(RESFUNC[,i],width = Width,space = 0,col = "darkgreen",
              beside =TRUE,border=NA,xlab="Duration",ylab=paste("R",i,sep="_"))
      axis(side=1,at=c(0:TPT))
    }
  }
  return(H)
}

####Function to evaluate EST, EFT, LST and LFT times of activity of a project####
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





