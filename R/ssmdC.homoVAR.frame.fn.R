#' @title  calculate the critical value of SSMD between a positive control and
#'   a negative reference in a data frame for an HTS study
#' @description A function to calculate the critical value of SSMD between a positive control and
#'   a negative reference in a data frame for an HTS study
#' @param  dataIn.df: a data frame with plate ID (or name), well usage and intensity
#' @param   negREF: name of wells to be used as a negative reference
#' @param   positiveCTRL: name of wells to be used as a positive control
#' @param   pREFtrim: the total portion of data to be trimmed in the negative reference
#' @param   Alpha: specify the significant level
#' @param   Beta: the population value of SSMD
#'
#' @author Xiaohua Douglas Zhang and Dandan Wang        07/2020
#' @return a matrix containing estimated SSMD value and variance for each plate
#' @examples  data("data.CVB3CRISPR", package="qcSSMDhomo")
#'   data.df = dataCVB3.df
#'   ssmdC.homoVAR.frame.fn(dataIn.df=data.df, negREF="Sample", positiveCTRL="Positive Control",
#'                          pREFtrim=0.05, Alpha=0.05, Beta=3)
#'
#' @export


ssmdC.homoVAR.frame.fn =
  function(dataIn.df, negREF="Sample", positiveCTRL="positive control", pREFtrim=0, Alpha=0.05, Beta=3)
  {
    #*****************************************************************************
    # function to calculate the critical value of SSMD between a positive control and
    #   a negative reference in a data frame for an HTS study
    # Author: Xiaohua Douglas Zhang and Dandan Wang        07/2020
    # Input
    #   dataIn.df: a data frame with plate ID (or name), well usage and intensity
    #   negREF: name of wells to be used as a negative reference
    #   positiveCTRL: name of wells to be used as a positive control
    #   pREFtrim: the total portion of data to be trimmed in the negative reference
    #   Alpha: specify the significant level
    #   Beta: the population value of SSMD
    # Output
    #   a matrix containing estimated SSMD value and variance for each plate.
    # Example:
	#   library(qcSSMD)
	#   data("data.CVB3CRISPR", package="qcSSMD")
    #   data.df = dataCVB3.df
	#   ssmdC.homoVAR.frame.fn(dataIn.df=data.df, negREF="Sample", positiveCTRL="Positive Control",
	#                          pREFtrim=0.05, Alpha=0.05, Beta=3)
    #*****************************************************************************

    plateID.vec = dataIn.df[,1]
    plateUniq.vec = unique(plateID.vec)
    nPlate = length( plateUniq.vec )
    N = length(plateID.vec)

    ssmdC.mat = matrix(NA, nrow=nPlate, ncol=3)
    for( i in 1:nPlate) {
      theData.df = dataIn.df[ plateID.vec == plateUniq.vec[i], ]
      inten.vec = theData.df[, 3]
      wellusage = theData.df[, 2]
      x = inten.vec[wellusage == negREF]
      y = inten.vec[wellusage == positiveCTRL]

      x = x[ !is.na(x) ]
      xSort = sort(x)
      nS1 = length(x)
      nRm = round(nS1*pREFtrim/2)
      xTrim.vec = xSort[ (nRm+1): (nS1-nRm) ]
      n1 = nS1 - nRm*2
      y = y[ !is.na(y) ]
      n2 = length(y)

      ssmdC.mat[i, ] = c(n1, n2, ssmdC.homoVAR.fn(n1, n2, Alpha, Beta))
    }
    dimnames(ssmdC.mat)= list(plateUniq.vec, c("nNeg", "nPos", "ssmdC"))
    return( ssmdC.mat )
  }
