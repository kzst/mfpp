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
tpq<- function(DSM,PEM, q)
{
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop(
      "Package \"pracma\" must be installed to use this function.",
      call. = FALSE
    )
  }
  TPQ=0  # Total Project Quality
  TPS=0  # Total Project Score (additive scores are assumed)

  for (i in 1:ncol(DSM))
  {
    if (DSM[i,i]>0)
    {TPS=TPS+PEM[i,i]}
    if (TPQ==0)
    {TPQ=1}

    TPQ=TPQ*(q[i]^PEM[i,i])
  }
  if (TPS>0)
    TPQ=(maxscore_PEM(DSM,PEM,pracma::ones(nrow(PEM))-PEM))*(TPQ^{1/TPS})
  return(TPQ)
}
