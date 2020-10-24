#' @title calculate UMVUE estimate of SSMD between two groups under homoscedasticity
#' @description A function to calculate UMVUE estimate of SSMD between two groups under homoscedasticity
#' @param  x the values of response in the first group (i.e., the reference group)
#' @param  y the values of response in the second group
#' @param  approx approx=TRUE, indicating the use of approximation for the coefficient in the UMVUE estimate of SSMD
#' @param  method specify the method for estimating SSMD. It must be either 'UMVUE', 'MLE' or'MM'.
#' @author Xiaohua Douglas Zhang originally in 2007, revised in 2020
#' @return a vector
#' @examples x = rnorm(320, mean=0, sd=1)
#'   y = rnorm(8, mean=3*sqrt(2), sd=1)
#'   ssmdEst.homoVAR.fn( x, y, pXtrim=0.05, approx=F, method="UMVUE")
#' @importFrom stats sd
#' @export





ssmdEst.homoVAR.fn = function( x, y, pXtrim=0, approx=FALSE, method="UMVUE")
{
  #*****************************************************************************************
  # A function to calculate UMVUE estimate of SSMD between two groups under homoscedasticity
  # Author: Xiaohua Douglas Zhang, originally in 2007, revised in 2020
  # Augment
  #   x: the values of response in the first group (i.e., the reference group)
  #   y: the values of response in the second group
  #   approx: approx=TRUE, indicating the use of approximation for the coefficient
  #           in the UMVUE estimate of SSMD
  #   method: specify the method for estimating SSMD. It must be
  #           either 'UMVUE', 'MLE' or'MM'. "
  # Example:
  #   x = rnorm(320, mean=0, sd=1)
  #   y = rnorm(8, mean=3*sqrt(2), sd=1)
  #   ssmdEst.homoVAR.fn( x, y, pXtrim=0.05, approx=FALSE, method="UMVUE")
  #*****************************************************************************************
  x = x[ !is.na(x) ]
  xSort = sort(x)
  nS1 = length(x)
  nRm = round(nS1*pXtrim/2)
  xTrim.vec = xSort[ (nRm+1): (nS1-nRm) ]
  n1 = nS1 - nRm*2
  m1 = mean(xTrim.vec)
  sd1 = sd(xTrim.vec)
  y = y[ !is.na(y) ]
  n2 = length(y)
  m2 = mean(y)
  sd2 = sd(y)
  ssmdEst.homoVAR.vec = ssmdEstCore.homoVAR.fn( m1, m2, sd1, sd2, n1, n2, approx=approx, method=method)
  return( ssmdEst.homoVAR.vec )
}
