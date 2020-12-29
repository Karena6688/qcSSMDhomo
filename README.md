
# qcSSMDhomo
#**************************************************************************************

#The main codes for using the R package qcSSMDhomo to generate the 

#tables, figures and values shown in the following paper.

#Zhang XD*#, Wang D#, Sun S, Zhang H*. Issues of z-factor and an 

#approach to avoid them for quality control in high-throughput screening studies 

#**************************************************************************************

install.packages("devtools") # If not already installed

devtools::install_github("Karena6688/qcSSMDhomo")

library(qcSSMDhomo)


#########################################################################################

## For Table 2 based on supplementary table 1 and supplementary figure 1

#########################################################################################

##*********************************************************************

## Supplementary Table 1: Critical Value

n1.vec = c(304, 12, 76, 76, 1276, 16,  8, 8, 4, 6, 5, 4, 4)

n2.vec = c(8,   8,   6,  3,   12, 12, 16, 8, 4, 2, 2, 3, 2)

ssmdC.vec = rep( NA, length(n1.vec) )

for( i in 1:length(n1.vec) ) {

  n1 = n1.vec[i]
  
  n2 = n2.vec[i]
  
  ssmdC.vec[i] = ssmdC.homoVAR.fn(n1, n2, Alpha=0.05, Beta=-3)
  
}

rbind(n1.vec, n2.vec, round(ssmdC.vec, 2) )

#n1.vec 304.00 12.00 76.00 76.00 1276.00 20.00  8.00  8.00  4.00  6.00  5.00   4.0  4.00

#n2.vec   8.00  8.00  6.00  3.00   12.00 12.00 16.00  8.00  4.00  2.00  2.00   3.0  2.00

#-       -3.47 -4.14 -3.66 -3.83   -3.35 -3.92 -4.03 -4.31 -5.21 -5.28 -5.54  -5.5 -5.93

##*********************************************************************

## For Supplementary Figure 1

#a 

n1.vec=c(76,76,76,76,76, 304,304,304,304,304,304,304,304,304, 1276,1276,1276,1276,1276,1276,1276,1276,1276,1276,1276)

n2.vec=c( 2, 3, 4, 5, 6,   4,  5,  6,  7,  8,  9, 10, 11,12,     4,   5,   6,   7,   8,  9,   10,  12,  14,  16, 20)

nT1=5; nT2=9; nT3=11

ssmdC.vec = rep( NA, length(n1.vec) )

for( i in 1:length(n1.vec) ) {

  n1 = n1.vec[i]

  n2 = n2.vec[i]

  ssmdC.vec[i] = ssmdC.homoVAR.fn(n1, n2, Alpha=0.05, Beta=3)

}

rbind(n1.vec, n2.vec, round(ssmdC.vec, 2) )

par(mar=c(7.6, 4.1, 1.1, 1.1))

caseLab.vec = paste( "n2=", n2.vec, ", n1=", n1.vec, sep="")

Nn = length(ssmdC.vec)

plot( c(1, Nn), range(c(3, 6.5, ssmdC.vec)), type="n", 
      xlab="", ylab="SSMD critical value", axes=FALSE)

axis(2, las=2); box()

axis(1, at=1:Nn, labels=caseLab.vec, las=2)

lines( 1:nT1, ssmdC.vec[1:nT1], col=1 )

lines( 1:nT1, rep(3.8, nT1), col=3 )

lines( 1:nT2+nT1, ssmdC.vec[1:nT2+nT1], col=1 )

lines( 1:nT2+nT1, rep(3.5, nT2), col=2 )

lines( 1:nT3+nT1+nT2, ssmdC.vec[1:nT3+nT2+nT1], col=1 )

lines( 1:nT3+nT2+nT1, rep(3.4, nT3), col=4 )

lines( c(-10, 1000), rep(3, 2), col="grey", lty=1, lwd=2 )

legend("topright", legend=c("QC cutoff 3.8 for the 96-well format", 
                            "QC cutoff 3.5 for the 384-well format", 
							"QC cutoff 3.4 for the 1536-well format"), 							
       col=c(3, 2, 4), lty=1 )

#b

#b          96                                304

n1.vec = c(3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6,  rep(8,6),  rep(16, 7),  rep(24, 7) )

n2.vec = c(2, 3, 2, 3, 4, 2, 3, 4, 2, 3, 4, 5,      3:8,        2:8*2,   2:8*2)

nT1=12; nT2=6+7; nT3=7

ssmdC.vec = rep( NA, length(n1.vec) )

for( i in 1:length(n1.vec) ) {

  n1 = n1.vec[i]

  n2 = n2.vec[i]

  ssmdC.vec[i] = ssmdC.homoVAR.fn(n1, n2, Alpha=0.05, Beta=3)

}

rbind(n1.vec, n2.vec, round(ssmdC.vec, 2) )

par(mar=c(7.6, 4.1, 1.1, 1.1))

caseLab.vec = paste( "n2=", n2.vec, ", n1=", n1.vec, sep="")

Nn = length(ssmdC.vec)

plot( c(1, Nn), range(c(3, 6, ssmdC.vec)), type="n", 
      xlab="", ylab="SSMD critical value", axes=FALSE)

axis(2, las=2); box()

axis(1, at=1:Nn, labels=caseLab.vec, las=2)

lines( c(-10, 1000), rep(3, 2), col="grey", lty=1, lwd=2 )

lines( 1:nT1, ssmdC.vec[1:nT1], col=1 )

lines( 1:nT1, rep(5.5, nT1), col=3 )

lines( 1:(nT2+nT3)+nT1, ssmdC.vec[1:(nT2+nT3)+nT1], col=1 )

lines( 1:nT2+nT1, rep(4.2, nT2), col=2 )

lines( 19:(nT1+nT2+nT3), rep(3.9, length(19:(nT1+nT2+nT3)) ), col=4 )

legend("topright", legend=c("QC cutoff 5.5 for the 96-well format", 
                            "QC cutoff 4.2 for the 384-well format", 
							"QC cutoff 3.9 for the 1536-well format"), 
       col=c(3, 2, 4), lty=1 )

########################################################################################
## Figure 2 (CVB3 CRISPR/CAS9 primary screen) Figure 2d 2e 2f
########################################################################################

#dataCVB3.df = read.csv("data.CVB3.csv")

#save(dataCVB3.df, file="data.CVB3CRISPR.RData")

data("data.CVB3CRISPR", package="qcSSMDhomo")

data.df = dataCVB3.df

data.df[1:3, ]

#-   plate wellType row column intensity

#1 plate01     MISC   1      1        NA

#2 plate01     MISC   1      2        NA

#3 plate01   Sample   1      3 -1.123124

plateName = "plate"; rowName = "row"; colName = "column"

wellName = "wellType"

unique(data.df[, wellName])

#[1] MISC             Sample           Negative Control Positive Control

#Levels: MISC Negative Control Positive Control Sample

condt = data.df[, plateName] == "plate01"

table( data.df[condt, wellName] )

#-            MISC Negative Control Positive Control           Sample 

#-               9                5                2               80 

negName="Negative Control"; sampleName = "Sample"

strongPos = "Positive Control"

Intensity = "intensity" 

plate.vec = as.vector(data.df[, plateName])

plates = unique(plate.vec)

nPlate = length(plates) #[1] 20

nRow=8; nCol=12; nWell = nRow * nCol

uniqWells = unique(data.df[, wellName])

uniqWells

#### test homoscedasticity

#Bartlett’s test: It compares the variances of two or more groups. The data must be normally distributed in this test also.

#Fligner-Killeen test: A non-parametric test that is very robust against departures from normality.

pval.bartlett.negCTRL.vec = pval.fligner.negCTRL.vec = ratioSD.negCTRL.vec = rep(NA, nPlate)

for(i in 1:nPlate) {

  condt = plate.vec==plates[i]

  theData.df0 = data.df[condt, ]

  theWelltype = theData.df0[, "wellType"]

  condt2 = is.element(theWelltype, c(strongPos, negName))

  theData.df = theData.df0[condt2,]

  pval.bartlett.negCTRL.vec[i]=bartlett.test(intensity ~ wellType, data=theData.df)$p.value

  a = tapply( theData.df[,"intensity"], INDEX=theData.df[,"wellType"], sd)

  ratioSD.negCTRL.vec[i] = a[strongPos]/a[negName]

  pval.fligner.negCTRL.vec[i]=fligner.test(intensity ~ wellType, data=theData.df)$p.value

} 

 quantile(ratioSD.negCTRL.vec, probs=c(0, 0.025, 0.05, 1:9/10, 0.95, 0.975, 1))

#-      0%      2.5%        5%       10%       20%       30%       40%       50% 

#0.6000091 0.7569881 0.9139670 1.0385482 1.1549859 1.3759802 1.6523345 1.9155355 

#-     60%       70%       80%       90%       95%     97.5%      100% 

#2.2665002 2.8967713 3.6489016 4.5234806 6.5387113 7.3441055 8.1494997 

pCut = 0.05

c(sum( pval.fligner.negCTRL.vec <= pCut), sum( pval.bartlett.negCTRL.vec <= pCut) )

#[1] 7 2

pCut = 0.01

c(sum( pval.fligner.negCTRL.vec <= pCut), sum( pval.bartlett.negCTRL.vec <= pCut) )

#[1] 0 1

#####calculating SSMD estimated and critical values and z-factor

ssmdEst.mat = ssmdEst.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                       negREF=negName, positiveCTRL=strongPos, 
                                       pREFtrim=0, approx=FALSE, method="UMVUE")

ssmdC.mat = ssmdC.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                   negREF=negName, positiveCTRL=strongPos, 
                                   pREFtrim=0, Alpha=0.05, Beta=-3)

cbind(ssmdC.mat, ssmdEst.mat)

sum( ssmdC.mat[,"ssmdC"] < ssmdEst.mat[,"SSMDest"] ) #[1] 6

which(ssmdC.mat[,"ssmdC"] < ssmdEst.mat[,"SSMDest"])

#plate03 plate07 plate09 plate13 plate14 plate15 

#-     3       7       9      13      14      15 

which(ssmdEst.mat[,"SSMDest"]> -5.5)

#plate03 plate07 plate09 plate13 plate14 plate15 

#-     3       7       9      13      14      15

zFactor.mat = 
zFactor.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                 negREF=negName, positiveCTRL=strongPos, 
                 pREFtrim=0.05, is.homoVAR = TRUE)

sum(zFactor.mat[, "zFactor"] < 0) # 1

which(zFactor.mat[, "zFactor"] < 0)

#plate07 

#-     7

sum(zFactor.mat[, "zFactor"] < 0.5) # 8

which(zFactor.mat[, "zFactor"] < 0.5)

#plate03 plate07 plate09 plate11 plate13 plate14 plate15 plate17 

#-     3       7       9      11      13      14      15      17 

#Figure 2: Estimated SSMD and critical value (CVB3)

#Fig.2d0

par(mar=c(5.1, 4.4, 1.1, 1.1))

par(mfrow=c(1,1))

uniqWells

#[1] empty        Sample       negativeCTRL positiveCTRL

col.vec=c("black", "grey", "green","red")

wellplotOrders = 1:4

pch.vec = rep(1, length(col.vec))

x.df=cbind(plateWelltoX.fn(data.df[,c(plateName, rowName, colName)], nRow, nCol, byRow=FALSE),
            "wellType"=data.df[, wellName])

Intensity = "intensity"  

y.vec = data.df[, Intensity]

yrange =  range(y.vec, na.rm=T) 

plot( range(x.df$x), yrange, type="n", xlab="Plate Number (Plate-well series by column)",
      ylab="Normalized Intensity", axes=FALSE, main="")

axis(2, las=2)

axis(1, at=unique(x.df[,"plateOrder"]-1)*nWell, label=unique(x.df[,"plateOrder"]), las=1 )

box()

for( i in wellplotOrders ) {

  condt = x.df[,"wellType"] == uniqWells[i]

  points( x.df[condt,"x"], y.vec[condt], col=col.vec[i], cex=0.5)

}

legend("bottomright", legend=c("test sgRNAs", "Negative control", "Positive control"), 
       col=col.vec[c(2,3,4)], pch=1, cex=1)

#figure 2d

par(mar=c(5.1, 4.4, 1.1, 1.1))	

xat = 1:nPlate; #yat = 0:(-15)*2

yRange = range( c( ratioSD.negCTRL.vec, 1 ) ); yRange

plot( c(1, length(plates)), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="Ratio of standard deviation",  main="")

axis(1, at=xat, labels=xat, las=2); axis(2, las=2); box()

points( 1:length(plates), ratioSD.negCTRL.vec, col="black", cex=1.5)

lines( 1:length(plates), ratioSD.negCTRL.vec, col="black", cex=0.1 ) 

pLab.negCTRL.vec = ifelse(pval.bartlett.negCTRL.vec <= 0.01, "+", 
                          ifelse(pval.bartlett.negCTRL.vec <= 0.05, "*", "") )

text( 1:length(plates), ratioSD.negCTRL.vec, labels=pLab.negCTRL.vec, col="red", cex=1.5)

lines( c(1-10, length(plates)+10 ), rep(1, 2),  col="grey" ) 

legend("topright", legend=c("Estimated ratio","0.01 < p <= 0.05", "p <= 0.01"), 
       col=c("black", "red", "red"), pch=c("o",  "*", "+"), cex=1)

#Fig.2e

par(mar=c(5.1, 4.4, 1.1, 1.1))	

xat = c(1:20); yat = 0:(-15)*2

SSMD.mat = cbind("SSMDest"=ssmdEst.mat[, "SSMDest"], "ssmdC"=ssmdC.mat[, "ssmdC"])

yRange = c( min(SSMD.mat), 2 )

plot( c(1, length(plates)), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="SSMD", 
	  main="")

axis(1, at=xat, labels=xat, las=2)

axis(2, las=2); box()

points( 1:length(plates), ssmdEst.mat[, "SSMDest"], col="black" )

lines( 1:length(plates), ssmdEst.mat[, "SSMDest"], col="black", lty=1 ) 

lines( c(1-10, length(plates)+10 ), rep(-5.5, 2), col="red" ) 

legend("topright", legend="SSMD cutoff of -5.5", col="red", lty=1, pch=NA, cex=1)

#Fig.2f

col.vec1 = c("red", "black")

xat = c(1:20); yat = 2:(-8)/2

yRange = range(zFactor.mat[, "zFactor"])

plot( c(1, length(plates) ), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="z-factor", 
	  main="")

axis(1, at=xat, labels=xat, las=2)

axis(2, las=2); box()

points( 1:length(plates), zFactor.mat[, "zFactor"], col="black" )

lines( 1:length(plates), zFactor.mat[, "zFactor"], col="black" )  

lines( c(1-10, length(plates)+10 ), rep(0.5, 2), col="red" ) 

lines( c(1-10, length(plates)+10 ), rep(0, 2), col="lightblue" ) 

legend("bottomright", legend=c("Positive control", "z-factor = 0", "z-factor = 0.5"), 
       col=c("black", "lightblue", "red"), lty=1, pch=c(1,NA,NA), cex=1)


#########################################################################################

## For supplementary Figure 2 (Mucin primary screen)

#########################################################################################

#dataMucin.df = read.csv("data.mucinHTS.csv")

#save(dataMucin.df, file="data.mucinHTS.RData")

data("data.mucinHTS", package="qcSSMDhomo")

data.df = dataMucin.df

data.df[1:3, ]

#-   plate         wellType row column  intensity

#1 plate01  No Cell Control   1      1 -3.5306719

#2 plate01 Positive Control   1      2 -1.3234041

#3 plate01            MISC1   1      3  0.2775412

plateName = "plate"; rowName = "row"; colName = "column"

wellName = "wellType"

unique(data.df[, wellName])

#[1] "No Cell Control"  "Positive Control" "MISC1"            "Test siRNAs"     

#[5] "MISC2"            "MISC3"            "Negative Control" "MISC4"  

condt = data.df[, plateName] == "plate01"

table( data.df[condt, wellName] )

#-          MISC1            MISC2            MISC3            MISC4 

#-             16               16                8                8 

#Negative Control  No Cell Control Positive Control      Test siRNAs 

#-              8               16                8              304 

negName="Negative Control"; sampleName = "Test siRNAs"

strongPos = "No Cell Control"

weakerPos = "Positive Control"

Intensity = "intensity" 

plate.vec = as.vector(data.df[, plateName])

plates = unique(plate.vec)

nPlate = length(plates) #[1] 20

nRow=16; nCol=24; nWell = nRow * nCol

uniqWells = unique(data.df[, wellName])

uniqWells

pval.bartlett.sample.vec = pval.fligner.sample.vec = ratioSD.sample.vec = rep(NA, nPlate)

pval.bartlett.negCTRL.vec = pval.fligner.negCTRL.vec = ratioSD.negCTRL.vec = rep(NA, nPlate)

for(i in 1:nPlate) {

  condt = plate.vec==plates[i]

  theData.df0 = data.df[condt, ]

  theWelltype = theData.df0[, "wellType"]

  condt2 = is.element(theWelltype, c(strongPos, negName))

  theData.df = theData.df0[condt2,]

  pval.bartlett.negCTRL.vec[i]=bartlett.test(intensity ~ wellType, data=theData.df)$p.value

  a = tapply( theData.df[,"intensity"], INDEX=theData.df[,"wellType"], sd)

  ratioSD.negCTRL.vec[i] = a[strongPos]/a[negName]

  pval.fligner.negCTRL.vec[i]=fligner.test(intensity ~ wellType, data=theData.df)$p.value
  
  condt3.1 = is.element( theWelltype, strongPos )

  x.vec = theData.df0[, Intensity]

  cutoffs = quantile( x.vec[theWelltype == sampleName], probs=c(0.025, 0.975), na.rm=TRUE )

  condt3.2 = x.vec >= cutoffs[1] & x.vec <= cutoffs[2] & is.element(theWelltype, sampleName)

  condt3 = condt3.1 | condt3.2

  theData.df = theData.df0[condt3,]

  pval.bartlett.sample.vec[i]=bartlett.test(intensity ~ wellType, data=theData.df)$p.value

  a = tapply( theData.df[,"intensity"], INDEX=theData.df[,"wellType"], sd)

  ratioSD.sample.vec[i] = a[strongPos]/a[sampleName]

  pval.fligner.sample.vec[i]=fligner.test(intensity ~ wellType, data=theData.df)$p.value

} 

 quantile(ratioSD.negCTRL.vec, probs=c(0, 0.025, 0.05, 1:9/10, 0.95, 0.975, 1))

#-       0%       2.5%         5%        10%        20%        30%        40%        50% 

#0.05845157 0.08424076 0.09796906 0.11071123 0.12958202 0.18695739 0.22730642 0.24401459 

#-      60%        70%        80%        90%        95%      97.5%       100% 

#0.29846500 0.36976919 0.43230948 0.57810413 0.68114107 0.92043109 4.08526167 

pCut = 0.05

c(sum( pval.fligner.negCTRL.vec <= pCut), sum( pval.bartlett.negCTRL.vec <= pCut) )

#[1] 53 66

pCut = 0.01

c(sum( pval.fligner.negCTRL.vec <= pCut), sum( pval.bartlett.negCTRL.vec <= pCut) )

#[1] 3 12

 quantile(ratioSD.sample.vec, probs=c(0, 0.025, 0.05, 1:9/10, 0.95, 0.975, 1))

#-      0%      2.5%        5%       10%       20%       30%       40%       50% 

#0.5229661 0.5966793 0.8342320 0.8849042 1.0593345 1.2668767 1.4610724 1.6203533 

#-     60%       70%       80%       90%       95%     97.5%      100% 

#1.7457276 1.9185457 2.3897030 2.7520714 3.2095790 3.9397742 4.5296211 

pCut = 0.05

c(sum( pval.fligner.sample.vec <= pCut), sum( pval.bartlett.sample.vec <= pCut) )

#[1] 43 53

pCut = 0.01

c(sum( pval.fligner.sample.vec <= pCut), sum( pval.bartlett.sample.vec <= pCut) )

#[1] 35 45


#calculating SSMD estimated and critical values and z-factor

ssmdEstStrong.mat = ssmdEst.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                             negREF=negName, positiveCTRL=strongPos, 
                                             pREFtrim=0.05, approx=FALSE, method="UMVUE")

ssmdCstrong.mat = ssmdC.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                         negREF=negName, positiveCTRL=strongPos, 
                                         pREFtrim=0.05, Alpha=0.05, Beta=-3)

zFactorStrong.mat = zFactor.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                     negREF=negName, positiveCTRL=strongPos, pREFtrim=0.05, is.homoVAR=FALSE)								   

zFactorStrong.mat2 = zFactor.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                     negREF=negName, positiveCTRL=strongPos, pREFtrim=0.05, is.homoVAR=TRUE)								   
								   
ssmdEstWeaker.mat = ssmdEst.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                             negREF=negName, positiveCTRL=weakerPos, 
                                             pREFtrim=0.05, approx=FALSE, method="UMVUE")

ssmdCweaker.mat = ssmdC.homoVAR.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                         negREF=negName, positiveCTRL=weakerPos, 
                                         pREFtrim=0.05, Alpha=0.05, Beta=-3)

zFactorWeaker.mat = zFactor.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                     negREF=negName, positiveCTRL=weakerPos,  pREFtrim=0.05, is.homoVAR=TRUE)

zFactorWeaker.mat2 = zFactor.frame.fn(dataIn.df=data.df[, c("plate", "wellType", "intensity")], 
                                     negREF=negName, positiveCTRL=weakerPos,  pREFtrim=0.05, is.homoVAR=FALSE)
									 
									 
which(ssmdCstrong.mat[,"ssmdC"] < ssmdEstStrong.mat[,"SSMDest"])

#plate31 

#-    31

which(ssmdEstStrong.mat[,"SSMDest"]> -4.2)

#plate31 

#-    31

sum(ssmdCweaker.mat[,"ssmdC"] < ssmdEstWeaker.mat[,"SSMDest"]) #11

which(ssmdCweaker.mat[,"ssmdC"] < ssmdEstWeaker.mat[,"SSMDest"])

#plate01 plate14 plate31 plate53 plate58 plate67 plate68 plate69 plate72 plate74 plate79 

#-     1      14      31      53      58      67      68      69      72      74      79 

which(ssmdEstWeaker.mat[,"SSMDest"]> -4.2)

#plate01 plate14 plate31 plate53 plate58 plate67 plate68 plate69 plate72 plate74 plate79 

#-     1      14      31      53      58      67      68      69      72      74      79 

which(zFactorStrong.mat2[,"zFactor"]<0)

#plate31 

#-    31 

sum(zFactorStrong.mat[,"zFactor"]<0) #1

which(zFactorStrong.mat2[,"zFactor"]<0)

#plate31 

#-    31 

sum(zFactorStrong.mat2[,"zFactor"]<0.5) #10

which(zFactorStrong.mat2[,"zFactor"]<0.5)

#plate17 plate31 plate32 plate33 plate34 plate36 plate37 plate38 plate41 plate61 

#-    17      31      32      33      34      36      37      38      41      61 

sum(zFactorStrong.mat[,"zFactor"]<0.5) #25

which(zFactorStrong.mat[,"zFactor"]<0.5)

#plate01 plate05 plate09 plate15 plate17 plate19 plate21 plate22 plate30 plate31 

#-     1       5       9      15      17      19      21      22      30      31 

#plate32 plate33 plate34 plate35 plate36 plate37 plate38 plate40 plate41 plate42 

#-    32      33      34      35      36      37      38      40      41      42 

#plate43 plate46 plate53 plate61 plate69 

#-    43      46      53      61      69 

sum(zFactorWeaker.mat2[,"zFactor"]<0) # 7			 

which(zFactorWeaker.mat2[,"zFactor"]<0)

#plate14 plate31 plate67 plate69 plate72 plate74 plate79 

#-    14      31      67      69      72      74      79 

sum(zFactorWeaker.mat[,"zFactor"]<0) # 10			 

which(zFactorWeaker.mat[,"zFactor"]<0)

#plate14 plate31 plate53 plate58 plate67 plate68 plate69 plate72 plate74 plate79 

#-    14      31      53      58      67      68      69      72      74      79 

sum(zFactorWeaker.mat2[,"zFactor"]<0.5)  # 47

which(zFactorWeaker.mat2[,"zFactor"]<0.5)

#plate01 plate02 plate03 plate04 plate05 plate07 plate08 plate10 plate11 plate12 plate14 

#-     1       2       3       4       5       7       8      10      11      12      14 

#plate16 plate17 plate18 plate20 plate23 plate24 plate25 plate26 plate27 plate31 plate36 

#-    16      17      18      20      23      24      25      26      27      31      36 

#plate37 plate47 plate50 plate53 plate56 plate57 plate58 plate60 plate61 plate62 plate64 

#-    37      47      50      53      56      57      58      60      61      62      64 

#plate65 plate66 plate67 plate68 plate69 plate70 plate71 plate72 plate73 plate74 plate76 

#-    65      66      67      68      69      70      71      72      73      74      76 

#plate77 plate78 plate79 

#-    77      78      79 

sum(zFactorWeaker.mat[,"zFactor"]<0.5) #55

which(zFactorWeaker.mat[,"zFactor"]<0.5)

cbind(which(zFactorWeaker.mat[,"zFactor"]<0.5), which(zFactorWeaker.mat2[,"zFactor"]<0.5)) 

table(c(which(zFactorWeaker.mat[,"zFactor"]<0.5), which(zFactorWeaker.mat2[,"zFactor"]<0.5)) )

##1  2  3  4  5  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 31 36 37 

##2  2  2  2  2  2  2  1  2  2  2  1  2  1  2  2  2  1  2  1  2  2  2  2  2  2  2  2 

#47 48 49 50 53 56 57 58 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 76 77 78 79 

##2  1  1  2  2  2  2  2  2  2  2  1  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2 

table(table(c(which(zFactorWeaker.mat[,"zFactor"]<0.5), which(zFactorWeaker.mat2[,"zFactor"]<0.5)) ))

#1  2 

#8 47 
	 
#supplementary Figure 2: Estimated SSMD and critical value (Mucin)

#supplementary Fig. 2a0

par(mar=c(5.1, 4.4, 1.1, 1.1))

par(mfrow=c(1,1))

uniqWells

#[1] "No Cell Control"  "Positive Control" "MISC1"            "Test siRNAs"     

#[5] "MISC2"            "MISC3"            "Negative Control" "MISC4"  

col.vec=c("red", "blue", "lightblue", "grey", "white", "white", "green", "white")

wellplotOrders = c(4, 7, 1, 2)

pch.vec = rep(1, length(col.vec))

x.df=cbind(plateWelltoX.fn(data.df[,c(plateName, rowName, colName)], nRow, nCol, byRow=FALSE),
            "wellType"=data.df[, wellName])

Intensity = "intensity"  

y.vec = data.df[, Intensity]

yrange =  range(y.vec, na.rm=T) 

plot( range(x.df$x), yrange, type="n", xlab="Plate Number (Plate-well series by column)",
      ylab="Normalized Intensity", axes=FALSE, main="")

axis(2, las=2)

axis(1, at=unique(x.df[,"plateOrder"]-1)*nWell, label=unique(x.df[,"plateOrder"]), las=2 )

box()

for( i in wellplotOrders ) {

  condt = x.df[,"wellType"] == uniqWells[i]

  points( x.df[condt,"x"], y.vec[condt], col=col.vec[i], cex=0.5)

}

legend("bottomright", legend=c("No-cell control", "Positive control", "Test siRNAs", "Negative control"), 
       col=col.vec[c(1, 2, 4, 7)], pch=1, cex=0.8)

#supplementary fig. 2a

par(mar=c(5.1, 4.4, 1.1, 1.1))	

xat = 1:nPlate; #yat = 0:(-15)*2

yRange = range( c( ratioSD.negCTRL.vec, ratioSD.sample.vec, 1 ) ); yRange

plot( c(1, length(plates)), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="ratio of standard deviation",  main="")

axis(1, at=xat, labels=xat, las=2); axis(2, las=2); box()

points( 1:length(plates), ratioSD.negCTRL.vec, col="green", cex=1.5)

lines( 1:length(plates), ratioSD.negCTRL.vec, col="green", cex=0.1 ) 

points( 1:length(plates), ratioSD.sample.vec, col="black", cex=1.5 )

lines( 1:length(plates), ratioSD.sample.vec, col="black", cex=0.5) 

pLab.sample.vec = ifelse(pval.bartlett.sample.vec <= 0.01, "+", 
                          ifelse(pval.bartlett.sample.vec <= 0.05, "*", "") )

text( 1:length(plates), ratioSD.sample.vec, labels=pLab.sample.vec, col="red", cex=1.2)

pLab.negCTRL.vec = ifelse(pval.bartlett.negCTRL.vec <= 0.01, "+", 
                          ifelse(pval.bartlett.negCTRL.vec <= 0.05, "*", "") )

text( 1:length(plates), ratioSD.negCTRL.vec, labels=pLab.negCTRL.vec, col="red", cex=1.2)

lines( c(1-10, length(plates)+10 ), rep(1, 2),  col="grey" ) 

legend("top", legend=c("Sample wells", "Negative control","0.01 < p <= 0.05", "p <= 0.01"), 
       col=c("black", "green", "red", "red"), 
      lty=c(1, 1, NA, NA), pch=c("o", "o", "*", "+"), cex=1)

#supplementary fig. 2b

par(mar=c(5.1, 4.4, 1.1, 1.1))	

xat = c(1:79); yat = 0:(-15)*2

SSMD.mat = cbind("SSMDest"=ssmdEstStrong.mat[, "SSMDest"], "ssmdC"=ssmdCstrong.mat[, "ssmdC"])

yRange = c( min(SSMD.mat), 2 )

plot( c(1, length(plates)), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="SSMD", 
	  main="")

axis(1, at=xat, labels=xat, las=2); axis(2, at=yat, labels=yat, las=2); box()

points( 1:length(plates), ssmdEstStrong.mat[, "SSMDest"], col="red" )

lines( 1:length(plates), ssmdEstStrong.mat[, "SSMDest"], col="red" ) 

points( 1:length(plates), ssmdEstWeaker.mat[, "SSMDest"], col="blue" )

lines( 1:length(plates), ssmdEstWeaker.mat[, "SSMDest"], col="blue" ) 

lines( c(1-10, length(plates)+10 ), rep(-4.2, 2), col="black" ) 

legend("topright", legend=c("SSMD estimated value of No-cell Control", "SSMD estimated value of Positive control"), 
       col=c("red", "blue"), lty=1, pch=1, cex=1)

#supplementary Fig. 2c

col.vec1 = c("red", "black")

xat = c(1:79); yat = 2:(-8)/2

yRange = range(c(zFactorStrong.mat[, "zFactor"], zFactorWeaker.mat[, "zFactor"]))

plot( c(1, length(plates) ), yRange, type="n", axes=FALSE, 
      xlab="Plate Number", ylab="z-factor", 
	  main="")

axis(1, at=xat, labels=xat, las=2); axis(2, at=yat, labels=yat, las=2); box()

points( 1:length(plates), zFactorStrong.mat[, "zFactor"], col="red" )

lines( 1:length(plates), zFactorStrong.mat[, "zFactor"], col="red" )  

points( 1:length(plates), zFactorWeaker.mat[, "zFactor"], col="blue" )

lines( 1:length(plates), zFactorWeaker.mat[, "zFactor"], col="blue" )  

lines( c(1-10, length(plates)+10 ), rep(0.5, 2), col="lightgreen" ) 

lines( c(1-10, length(plates)+10 ), rep(0, 2), col="black" ) 

legend("bottomright", legend=c("No-cell Control", "Positive control", "z-factor = 0", "z-factor = 0.5"), 
       col=c("red", "blue", "black", "lightgreen"), lty=1, pch=c(1,1,NA,NA), cex=1)

###########################################################
####end of codes for the paper ######################


## Data

The related dataset is stored in Excel files in the subdirectory "extdata".

The data may also be obtained through running the following command lines in R after installing qcSSMD

library(qcSSMDhomo)

data("data.CVB3CRISPR", package="qcSSMDhomo")

The CVB3 CRISPR/CAS9 screen was published originally in

Kim, H. S. et al. Arrayed CRISPR screen with image-based assay reliably uncovers host genes required for coxsackievirus infection. Genome Res 28, 859-868 (2018).

## Functions

qcSSMDhomo contains the following major functions.

plateWelltoX.fn

ssmdC.homoVAR.fn

ssmdC.homoVAR.frame.fn

ssmdEst.homoVAR.fn

ssmdEst.homoVAR.frame.fn

ssmdEstCore.homoVAR.fn

zFactor.frame.fn

## Author for codes
Dandan Wang and Xiaohua Douglas Zhang, Ph.D., Professor, University of Macau

## Author for the Paper
Zhang XD*#, Wang D#, Sun S, Zhang H*. Issues of z-factor and an approach to avoid them for quality control in high-throughput screening studies. Bioinformatics (in press) 

