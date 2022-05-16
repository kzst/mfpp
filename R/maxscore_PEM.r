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
