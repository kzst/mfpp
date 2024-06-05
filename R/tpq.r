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
