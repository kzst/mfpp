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
  TPC<-0
  for (i in 1:length(CD)){
    TPC=TPC+as.numeric(DSM[i,i])*as.numeric(CD[i])
  }
  return(as.numeric(TPC))
}
