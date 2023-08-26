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
tpc <- function(DSM,CD){
  TPC<-sum(ceiling(as.numeric(diag(DSM)))*as.numeric(CD))
  return(as.numeric(TPC))
}
