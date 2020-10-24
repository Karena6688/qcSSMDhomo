
#' @title  calculate the plug-in value of z-factor between a positive control and
#   a negative reference in a data frame for an HTS study
#' @description A function to calculate the plug-in value of z-factor between a positive control and
#   a negative reference in a data frame for an HTS study
#' @param  dataIn.df a data frame with plate ID (or name), well usage and intensity
#' @param  negREF name of wells to be used as a negative reference
#' @param  positiveCTRL name of wells to be used as a positive control
#' @param  pREFtrim the total portion of data to be trimmed in the negative reference
#' @param  is.homoVAR indicator for whether to use homoscedasticity or heteroscedasticity
#' @param  k: the coefficient in z-factor
#'
#' @author Xiaohua Douglas Zhang and Dandan Wang        07/2020
#' @return a matrix containing estimated plug-in z-factor for each plate
#' @examples  data("data.CVB3CRISPR", package="qcSSMDhomo")
#'   data.df = dataCVB3.df
#'   zFactor.frame.fn(dataIn.df=data.df, negREF="Sample", positiveCTRL="Positive Control",
#'                    pREFtrim=0.05, pREFtrim=0, is.homoVAR=TRUE, k=3)
#'
#' @importFrom stats sd
#' @export


zFactor.frame.fn =
  function(dataIn.df, negREF="Sample", positiveCTRL="posCTRL", pREFtrim=0, is.homoVAR=TRUE, k=3)
  {
    #*****************************************************************************
    # function to calculate the plug-in value of z-factor between a positive control and
    #   a negative reference in a data frame for an HTS study
    # Author: Xiaohua Douglas Zhang and Dandan Wang        07/2020
    # Input
    #   dataIn.df: a data frame with plate ID (or name), well usage and intensity
    #   negREF: name of wells to be used as a negative reference
    #   positiveCTRL: name of wells to be used as a positive control
    #   pREFtrim: the total portion of data to be trimmed in the negative reference
    #   is.homoVAR: indicator for whether to use homoscedasticity or heteroscedasticity
    #   k: the coefficient in z-factor
    # Output
    #   a matrix containing estimated plug-in z-factor for each plate.
    # Example:
    #
    #*****************************************************************************

    plateID.vec = dataIn.df[,1]
    plateUniq.vec = unique(plateID.vec)
    nPlate = length( plateUniq.vec )
    N = length(plateID.vec)

    zFactor.mat = matrix(NA, nrow=nPlate, ncol=3)
    for( i in 1:nPlate) {
      theData.df = dataIn.df[ plateID.vec == plateUniq.vec[i], ]
      inten.vec = theData.df[, 3]
      wellusage = theData.df[, 2]

      x = inten.vec[wellusage == negREF]
      x = x[ !is.na(x) ]
      xSort = sort(x)
      nS1 = length(x)
      nRm = round(nS1*pREFtrim/2)
      xTrim.vec = xSort[ (nRm+1): (nS1-nRm) ]
      n1 = nS1 - nRm*2
      m1 = mean(xTrim.vec)
      sd1 = sd(xTrim.vec)

      y = inten.vec[wellusage == positiveCTRL]
      y = y[ !is.na(y) ]
      n2 = length(y)
      m2 = mean(y)
      sd2 = sd(y)

	  SDcombined = ifelse(is.homoVAR,
	                      2*sqrt( ((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2) ), sd1+sd2)
      zFactor.mat[i, ] = c( n1, n2, 1-k*SDcombined/abs(m1-m2) )
    }
    dimnames(zFactor.mat)= list(plateUniq.vec, c("nNeg", "nPos", "zFactor"))
    return( zFactor.mat )
  }
