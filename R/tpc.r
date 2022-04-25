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
tpc = function (DSM, C)
{
  TPC<-t(diag(ceiling(DSM)))%*%C
  rownames(TPC)<-"TPC"
  colnames(TPC)<-"TPC"
  return(TPC)
}
