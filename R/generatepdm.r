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
#' @importFrom stats na.omit
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
