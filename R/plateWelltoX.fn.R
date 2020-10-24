#' @title Convert plate, row and column numbers to a series of values
#' @description a function to generate value in x-axis based on plate number and well positions in a plate.
#' @param data.df consists of three columns: plate names or index, rows and columns
#' @param nRow number of rows in a plate
#' @param nCol number of columns in a plate
#' @param byRow whether the well series is arranged by row
#' @author Xiaohua Douglas Zhang   2017
#' @return a data frame
#' @examples plateID = rep(1:5, each=384)
#'   rowCore = rep(1:16, each=24); rows = rep(rowCore, 5)
#'   columnCore = rep(1:24, 16); columns = rep(columnCore, 5)
#'   data.df = data.frame(plateID, rows, columns)
#'   x.df = plateWelltoX.fn(data.df, nRow=16, nCol=24, byRow=TRUE)
#'
#' @export

plateWelltoX.fn = function(data.df, nRow=16, nCol=24, byRow=TRUE) {
#****************************************************************************
# function to generate value in x-axis based on plate number and
#   well positions in a plate.
# data.df: consists of three columns: plate names or index, rows and columns
# nRow: number of rows in a plate
# nCol: number of columns in a plate
# byRow: whether the well series is arranged by row
# Author: Xiaohua Douglas Zhang, 2017
# Example:
#   plateID = rep(1:5, each=384)
#   rowCore = rep(1:16, each=24); rows = rep(rowCore, 5)
#   columnCore = rep(1:24, 16); columns = rep(columnCore, 5)
#   data.df = data.frame(plateID, rows, columns)
#   x.df = plateWelltoX.fn(data.df, nRow=16, nCol=24, byRow=TRUE)
#****************************************************************************
  plate.vec = data.df[,1]
  platesUniq = unique(plate.vec)
  plateOrder.vec = rep( NA, length(plate.vec) )

  for( i in 1:length(platesUniq) ) plateOrder.vec[plate.vec==platesUniq[i]] = i
  if( byRow ) {
    x.vec = (plateOrder.vec-1)*nRow*nCol + (data.df[,2]-1)*nCol + data.df[,3]
  } else {
    x.vec = (plateOrder.vec-1)*nRow*nCol + (data.df[,3]-1)*nRow + data.df[,2]
  }
  return( data.frame(x=x.vec, plateOrder = plateOrder.vec) )
}

