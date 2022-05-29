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
plot.PDM_matrix <- function(x,...){
  if ("PDM_matrix" %in% class(x)){
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop(
        "Package \"igraph\" must be installed to use this function.",
        call. = FALSE
      )
    }
    PDM<-x
    class(PDM)<-"PDM_matrix"
    N<-dim(PDM)[1]
    M<-dim(PDM)[2]
    if (N>M){
      stop(
        "number of rows must be less or equal than the rows",
        call. = FALSE
      )
    }else{
      if (is.flexible(PDM)){
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
        diag(minpdm)<-0
        n<-dim(minpdm)[1]
        plot(igraph::graph.adjacency(minpdm[1:n,1:n]),main="Minimal Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        maxpdm<-truncpdm(maxPDM)
        diag(maxpdm)<-0
        n<-dim(maxpdm)[1]
        plot(igraph::graph.adjacency(maxpdm[1:n,1:n]),main="Maximal Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        minimaxpdm<-truncpdm(minimaxPDM)
        diag(minimaxpdm)<-0
        n<-dim(minimaxpdm)[1]
        plot(igraph::graph.adjacency(minimaxpdm[1:n,1:n]),main="Minimax Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        maximinpdm<-truncpdm(maximinPDM)
        diag(maximinpdm)<-0
        n<-dim(maximinpdm)[1]
        plot(igraph::graph.adjacency(maximinpdm[1:n,1:n]),main="Maximin Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

        mostpdm<-truncpdm(mostPDM)
        diag(mostpdm)<-0
        n<-dim(mostpdm)[1]
        plot(igraph::graph.adjacency(mostpdm[1:n,1:n]),main="Most-likely/Most-desired Structure",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")

    }else{
        pdm<-truncpdm(PDM)
        diag(pdm)<-0
        n<-dim(pdm)[1]
        plot(igraph::graph.adjacency(pdm[1:n,1:n]),main="Logic Network",
                     layout=igraph::layout_as_tree,vertex.shape="square",vertex.color="green")
      }

    }
  }else{
    plot(x,...)
  }
}
