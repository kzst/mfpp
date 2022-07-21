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
summary.PDM_const <- function(object, digits =  getOption("digits"), ...) {
  if ("PDM_const" %in% class(object)){
    cat("\nSummary of the PDM constraints structure:\n")
    if (!is.null(object$Ct)) cat("\nTime constraint (Ct): ",
                                 round(object$Ct,digits = digits))
    if (!is.null(object$Cc)) cat("\nConst constraint (Cc): ",
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
  if ("PDM_matrix" %in% class(object)){
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
  if ("PDM_list" %in% class(object)){
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
  if ("Set_PDM_matrix" %in% class(object)){
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
  if ("Set_PDM_list" %in% class(object)){
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
  if ("Collection_PDM" %in% class(object)){
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
  if ("TPT" %in% class(object)){
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
