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
plot.PDM_matrix <- function(x,w=getOption("w"),Rs=getOption("Rs"),...){
  if ("PDM_matrix" %in% class(x)){
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
    PDM<-x
    class(PDM)<-"PDM_matrix"
    N<-pracma::size(PDM,1)
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
        plot(g,main="Original Logic Network",
             layout=igraph::layout_as_tree,vertex.shape="square")

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


        if (!is.null(c)){
          plot(g,main="Minimal Structure",layout=igraph::layout_as_tree,
               vertex.shape="square",vertex.label=paste("d",igraph::V(g)$weight,sep="="))

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
               vertex.shape="square")

        }

        maxpdm<-truncpdm(maxPDM)
        n<-pracma::size(maxpdm,1)
        m<-pracma::size(maxpdm,2)
        c<-NULL
        if (!is.null(w)) { # Number of completion mode is specified
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


        if (!is.null(c)){
          plot(g,main="Maximal Structure",layout=igraph::layout_as_tree,
               vertex.shape="square",vertex.label=paste("d",
                                                igraph::V(g)$weight,sep="="))

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
               vertex.shape="square")

        }

        minimaxpdm<-truncpdm(minimaxPDM)
        diag(minimaxpdm)<-0
        n<-pracma::size(minimaxpdm,1)
        plot(igraph::graph.adjacency(minimaxpdm[1:n,1:n]),main="Minimax Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        maximinpdm<-truncpdm(maximinPDM)
        diag(maximinpdm)<-0
        n<-pracma::size(maximinpdm,1)
        plot(igraph::graph.adjacency(maximinpdm[1:n,1:n]),main="Maximin Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        mostpdm<-truncpdm(mostPDM)
        diag(mostpdm)<-0
        n<-pracma::size(mostpdm,1)
        plot(igraph::graph.adjacency(mostpdm[1:n,1:n]),main="Most-likely/Most-desired Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

    }else{
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
        plot(g,main="Logic Network",
                     layout=igraph::layout_as_tree,vertex.shape="square")

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
          pt.bg  = c("black", "red"),
          pch    = 45,
          cex    = 1,
          bty    = "n",
          title  = "Dependencies"
        )
      }

    }
  }else{
    plot(x,...)
  }
}
