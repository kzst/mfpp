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
    oldpar<-graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow=c(1,1))
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
          graphics::legend(
            "topleft",
            legend = c("mandatory", "supplementary"),
            pt.bg  = c("green", "yellow"),
            pch    = 22,
            cex    = 1,
            bty    = "n",
            title  = "Tasks"
          )
          graphics::legend(
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

            graphics::legend(
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

            graphics::legend(
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
        oldpar<-graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar))
        graphics::par(mfrow=c(1,n))
      }
      if (!is.null(minCONST$Ct)&&!is.null(maxCONST$Ct))
        graphics::barplot(cbind(minCONST$Ct,maxCONST$Ct),
                names.arg = c("min_C_t","max_C_t"),
                ylab = "TPT",main = "Duration constraints",
                col=col)
      if (!is.null(minCONST$Cc)&&!is.null(maxCONST$Cc))
        graphics::barplot(cbind(minCONST$Cc,maxCONST$Cc),
                names.arg = c("min_C_c","max_C_c"),
                ylab = "TPC",main = "Cost constraints",
                col=col)
      if (!is.null(minCONST$Cq)&&!is.null(maxCONST$Cq))
        graphics::barplot(cbind(minCONST$Cq,maxCONST$Cq),
                names.arg = c("min_C_q","max_C_q"),
                ylab = "TPQ",main = "Quality constraints",
                col=col)
      if (!is.null(minCONST$Cs)&&!is.null(maxCONST$Cs))
        graphics::barplot(cbind(minCONST$Cs,maxCONST$Cs),
                names.arg = c("min_C_s","max_C_s"),
                ylab = "TPS",main = "Scope/score constraints",
                col=col)
      if (!is.null(minCONST$CR)&&!is.null(maxCONST$CR))
        graphics::barplot(cbind(minCONST$CR,maxCONST$CR),
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
