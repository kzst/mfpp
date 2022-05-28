#-----------------------------------------------------------------------------#
#                                                                             #
#  RISK-BASED MULTIVARIATE CONTROL CHARTS                                     #
#                                                                             #
#  Written by: Aamir Saghir, Attila I. Katona, Zsolt T. Kosztyan              #
#              Department of Quantitative Methods                             #
#              University of Pannonia, Hungary                                #
#              kzst@gtk.uni-pannon.hu                                         #
#                                                                             #
# Last modified: March 2022                                                   #
#-----------------------------------------------------------------------------#

#' @export
summary.mfpp <- function(object, digits =  getOption("digits"), w=NULL,
                         Rs=NULL, ...) {

  if ("PDM_const" %in% class(object)){
    cat("\nSummary of the PDM constraints structure:\n")
    if (!is.null(object$Ct)) cat("\nTime constraint (Ct): ",
                                 round(object$Ct,digits = digits))
    if (!is.null(object$Cc)) cat("\nConst constraint (Cc): ",
                                 round(object$Cc,digits =digits))
    if (!is.null(object$Cs)) cat("\nScore/scope constraint (Cs): ",
                                 round(object$Cs,digits =digits))
    if (!is.null(object$Cq)) cat("\nQuailty constraint (Cq): ",
                                 round(object$Cq,digits =digits))
    if (!is.null(object$CR)) {
      cat("\nResource constraint(s) (CR):\n")
      round(object$CR,digits=digits)
    }
  }else{
    if ("PDM_matrix" %in% class(object)){
      cat("\nPrint PDM matrix:\n")
      print(object,digits=digits)
      if ((!is.null(w))||(!is.null(Rs))){
        maxCONST<-percent(object,type=c("c","q","qd","r","s","t"),w=w,Rs=Rs,
                          ratio=1)
        minCONST<-percent(object,type=c("c","q","qd","s","t"),w=w,Rs=Rs,
                          ratio=0)
        cat("\nMinimal constraints:\n")
        summary.mfpp(minCONST,digits=digits)
        cat("\n\nMaximal constraints:\n")
        summary.mfpp(maxCONST,digits=digits)
      }else{
        maxCONST<-percent(object,type=c("s"),ratio=1)
        minCONST<-percent(object,type=c("s"),ratio=0)
        cat("\nMinimal constraints:\n")
        summary.mfpp(minCONST,digits=digits)
        cat("\n\nMaximal constraints:\n")
        summary.mfpp(maxCONST,digits=digits)
        cat("\n\n")
        warning("Number of completion modes (w), and ",
        "number of resources (Rs) should be specified ",
        "to calculate constraints of demands.")
      }
    }else{
      if ("PDM_list" %in% class(object)){
        cat("\nPrint PDM list:\n")
        if (!is.null(object$w)) cat("\nNumber of completion modes (w): ",
                                    object$w)
        if (!is.null(object$Rs)) cat("\nNumber of resources (Rs): ",
                                    object$Rs)
        summary.mfpp(object$PDM,digits=digits,w=object$w,Rs=object$Rs)
      }else{
        if (("Set_PDM_matrix" %in% class(object))||
            ("Set_PDM_list" %in% class(object))){
          cat("\nSummary of main structures:\n")
          if (!is.null(object$minstruct)) {
            cat("\n\n\nSummary of minimal structure:\n")
            summary.mfpp(object$minstruct,digits = digits)
          }
          if (!is.null(object$maxstruct)) {
            cat("\n\n\nSummary of maximal structure:\n")
            summary.mfpp(object$maxstruct,digits = digits)
          }
          if (!is.null(object$minimaxstruct)) {
            cat("\n\n\nSummary of minimax structure:\n")
            summary.mfpp(object$minimaxstruct,digits = digits)
          }
          if (!is.null(object$maximinstruct)) {
            cat("\n\n\nSummary of maximin structure:\n")
            summary.mfpp(object$maximinstruct,digits = digits)
          }
          if (!is.null(object$moststruct)) {
            cat("\n\n\nSummary of most likely/most desired structure:\n")
            summary.mfpp(object$moststruct,digits = digits)
          }
        }else{
          summary(object,digits,...)
        }
      }
    }
  }
}
