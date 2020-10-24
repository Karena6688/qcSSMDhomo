#' @title calculate UMVUE estimate of SSMD between a positive control and
#'   a negative reference in a data frame for an HTS study
#' @description A function to calculate UMVUE estimate of SSMD between a positive control and
#'   a negative reference in a data frame for an HTS study
#' @param  dataIn.df a data frame with plate ID (or name), well usage and intensity
#' @param  negREF name of wells to be used as a negative reference
#' @param  positiveCTRL name of wells to be used as a positive control
#' @param  pREFtrim the total portion of data to be trimmed in the negative reference
#' @param  approx approx=TRUE, indicating the use of approximation for the coefficient
#'           in the UMVUE estimate of SSMD
#' @param   method specify the method for estimating SSMD. It must be
#'           either 'UMVUE', 'MLE' or'MM'. "
#' @author Xiaohua Douglas Zhang and Dandan Wang        07/2020
#' @return a matrix containing estimated SSMD value and variance for each plate
#' @examples  data("data.CVB3CRISPR", package="qcSSMDhomo")
#'   data.df = dataCVB3.df
#'   ssmdEst.homoVAR.frame.fn(dataIn.df=data.df, negREF="Sample", positiveCTRL="Positive Control",
#'                          pREFtrim=0.05, Alpha=0.05, Beta=3)
#'
#'
#' @export



ssmdEst.homoVAR.frame.fn =
  function(dataIn.df, negREF="Sample", positiveCTRL="positive control",
           pREFtrim=0, approx=FALSE, method="UMVUE")
  {
    #*****************************************************************************
    # function to calculate UMVUE estimate of SSMD between a positive control and
    #   a negative reference in a data frame for an HTS study
    # Author: Xiaohua Douglas Zhang and Dandan Wang        07/2020
    # Input
    #   dataIn.df: a data frame with plate ID (or name), well usage and intensity
    #   negREF: name of wells to be used as a negative reference
    #   positiveCTRL: name of wells to be used as a positive control
    #   pREFtrim: the total portion of data to be trimmed in the negative reference
    #   approx: approx=TRUE, indicating the use of approximation for the coefficient
    #           in the UMVUE estimate of SSMD
    #   method: specify the method for estimating SSMD. It must be
    #           either 'UMVUE', 'MLE' or'MM'. "
    # Output
    #   a matrix containing estimated SSMD value and variance for each plate.
    #*****************************************************************************
    plateID.vec = dataIn.df[,1]
    plateUniq.vec = unique(plateID.vec)
    nPlate = length( plateUniq.vec )
    N = length(plateID.vec)

    ssmdEst.mat = matrix(NA, nrow=nPlate, ncol=2)
    for( i in 1:nPlate) {
      theData.df = dataIn.df[ plateID.vec == plateUniq.vec[i], ]
      inten.vec = theData.df[, 3]
      wellusage = theData.df[, 2]
      negInten.vec = inten.vec[wellusage == negREF]
      posInten.vec = inten.vec[wellusage == positiveCTRL]
      ssmdEst.mat[i,] = ssmdEst.homoVAR.fn( negInten.vec, posInten.vec, pREFtrim, approx, method)
    }
    dimnames(ssmdEst.mat)= list(plateUniq.vec, c("SSMDest", "SSMDvar"))
    return( ssmdEst.mat )
  }
