#' @title calculate UMVUE estimate of SSMD
#' @description A function to calculate UMVUE estimate of SSMD
#' @param m1 sample mean of values in the first group (i.e., the reference group)
#' @param m2 sample mean of values in the second group
#' @param s1 sample standard deviation of values in the first group (i.e., the reference group)
#' @param s2 sample standard deviation of values in the second group
#' @param n1 sample size in the first group (i.e., the reference group)
#' @param n2 sample size in the second group
#' @param approx approx=TRUE, indicating the use of approximation for the coefficient in the UMVUE estimate of SSMD
#' @param method: specify the method for estimating SSMD. It must be 
#' @author Xiaohua Douglas Zhang  originally in 2007, revised in 2020
#' @return a vector
#' @examples  x = rnorm(304, mean=0, sd=1)
#'   y = rnorm(8, mean=3*sqrt(2), sd=1)
#'   ssmdEstCore.homoVAR.fn( mean(x), mean(y), sd(x), sd(y), length(x), length(y), method="UMVUE") 
#'
#' @export


## Function for calculate estimated value of SSMD for the general rule
## in Table 2
ssmdEstCore.homoVAR.fn = 
  function( m1, m2, s1, s2, n1, n2, approx=FALSE, method="UMVUE")
  {
    #*****************************************************************************
    # A function to calculate UMVUE estimate of SSMD
    # Author: Xiaohua Douglas Zhang, originally in 2007, revised in 2020
    # Augment
    #   m1: sample mean of values in the first group (i.e., the reference group)
    #   m2: sample mean of values in the second group
    #   s1: sample standard deviation of values in the first group (i.e., the reference group)
    #   s2: sample standard deviation of values in the second group
    #   n1: sample size in the first group (i.e., the reference group)
    #   n2: sample size in the second group
    #   approx: approx=TRUE, indicating the use of approximation for the coefficient
    #           in the UMVUE estimate of SSMD
    #   method: specify the method for estimating SSMD. It must be 
    #           either 'UMVUE', 'MLE' or'MM'. "
    # Example:
    #   x = rnorm(304, mean=0, sd=1)
    #   y = rnorm(8, mean=3*sqrt(2), sd=1)
    #   ssmdEstCore.homoVAR.fn( mean(x), mean(y), sd(x), sd(y), length(x), length(y), method="UMVUE") 
    #*****************************************************************************
    SS1 = ifelse( n1==1, 0, (n1-1)*s1^2 )
    SS2 = ifelse( n2==1, 0, (n2-1)*s2^2 )
    N = n1 + n2
    if( N < 4 ) { 
      stop("Error: there must be at least 4 non-missing values in total to estimate SSMD. ")	
    }
    if( approx ) { K = N - 3.5 } else {
      K = 2* (exp(lgamma( (N-2)/2 ) - lgamma( (N-3)/2 )))^2
    }
    varSSMD0 = N*K/(2*n1*n2*(N-4)) + (K-N+4)/(N-4)* N/2 * (m2-m1)^2/(SS1 + SS2)
    if( method == "UMVUE" ) {
      SSMDest = (m2 - m1)/sqrt( 2/K *(SS1 + SS2) )
      varSSMD = varSSMD0
    } else if( method == "MLE" ) {
      SSMDest = (m2 - m1)/sqrt( 2/N *(SS1 + SS2) )
      varSSMD = N/K * varSSMD0
    } else if( method == "MM" ) {
      SSMDest = (m2 - m1)/sqrt( 2/(N-2)*(SS1 + SS2) )
      varSSMD = (N-2)/K * varSSMD0
    } else {
      print("Error: 'method' must be either 'UMVUE', 'MLE' or'MM'. ")
    }
    ssmdEst.vec = c("SSMDest" = SSMDest, "SSMDvar" = varSSMD) 
    return( ssmdEst.vec)
  }