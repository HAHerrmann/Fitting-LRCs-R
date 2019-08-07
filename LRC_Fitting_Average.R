# Installing and loading the required packages
if (!require("Hmisc")) install.packages("Hmisc")
if (!require("readxl")) install.packages("readxl")
library(readxl)
library(Hmisc)

# Reading in the data
dat = read_excel("LightResponseCurvesEx.xlsx", sheet = 1)

# Averaging the LR curves 
LRcurves <- function(geno,temp,day){
  select = dat[which(dat$Genotype==geno & dat$Temperature==temp & dat$Day==day),]
  LRcurve = c()
  LRerror = c()
  for (i in c(0,125,250,500,750,1000,1500,2000)){
    x = which(select$Intensity==i)
    val = mean(select$A[x])
    err = sd(select$A[x])
    LRcurve = c(LRcurve,val)
    LRerror = c(LRerror,err)
  }
  return(list("LRcurve"=LRcurve,"LRerror"=LRerror))
}

# Predicting LR curve functions from the data 
LRcurve_predict <- function(P){
  curve.nlslrc = nls(P ~ (1/(2*theta))*(phi*I+Pmax-sqrt((phi*I+Pmax)^2-4*phi*theta*Pmax*I))-Rd,
                     start=list(Pmax=(max(P)-min(P)),phi=0.05,Rd=-min(P),theta=1)) 
  return(curve.nlslrc)
}

# Plotting and fitting Col0 Control0 as an example
P20=LRcurves("Col0",20,0)$LRcurve 
P20_E=LRcurves("Col0",20,0)$LRerror
I=c(0,125,250,500,750,1000,1500,2000)
P20_Pred=LRcurve_predict(P20) 
par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
print(plot(I,P20,xlab="", ylab="", ylim=c(-2,20),cex.lab=1.2,cex.axis=1.5,cex=2,pch=5))
errbar(I,P20,P20+P20_E,P20-P20_E,add=TRUE,col='black')
lines(0:2000,predict(P20_Pred, newdata=data.frame(I=0:2000)),lwd=1)
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
mtext(expression(P*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
title("Light Response Curve of Arabidopsis thaliana",line=-3,cex.main=2)
