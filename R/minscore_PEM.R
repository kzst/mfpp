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
