#' @title calculate the critical value of SSMD
#' @description A function to calculate the critical value of SSMD
#' @param n1 a value or a vector for the first group (i.e., the reference group)
#' @param n2 a vector for the second group
#' @param Alpha specify the significant level
#' @param Beta the population value of SSMD
#' @author Xiaohua Douglas Zhang   11/2019
#' @return a value
#' @examples ssmdC.homoVAR.fn(304, 16, Alpha=0.05, Beta=3)
#' @importFrom stats qt
#' @export




## Function for calculate critical value of SSMD for the general rule
## in Table 2
ssmdC.homoVAR.fn = function(n1, n2, Alpha=0.05, Beta=3)
  #*****************************************************************************
  # A function to calculate the critical value of SSMD
  # Author: Xiaohua Douglas Zhang        11/2019
  # Augment
  #   n1: a value or a vector for the first group (i.e., the reference group)
  #   n2: a vector for the second group
  #   Alpha: specify the significant level
  #   Beta: the population value of SSMD
  # Example:
  #   ssmdC.homoVAR.fn(304, 16, Alpha=0.05, Beta=3)
  #******************************************************************************
{
  DF = n1+n2-2
  NCPpart = sqrt( 2/( 1/n1 + 1/n2 ) )
  Tcritical = ifelse( Beta > 0, qt(p=1-Alpha, df=DF, ncp=Beta*NCPpart),
                      qt(p=Alpha, df=DF, ncp=Beta*NCPpart) )
  ssmdCritical = sqrt(2/DF)*exp(lgamma(DF/2)-lgamma((DF-1)/2))/NCPpart*Tcritical
  return(ssmdCritical)
}
