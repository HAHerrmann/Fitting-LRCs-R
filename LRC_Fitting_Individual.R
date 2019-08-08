

# Predicting LR curve functions from the data 
LRcurve_predict <- function(P,phi_guess=0.05,traceYN=FALSE){
  curve.nlslrc = nls(P ~ (1/(2*theta))*(phi*I+Pmax-sqrt((phi*I+Pmax)^2-4*phi*theta*Pmax*I))-Rd,
                     start=list(Pmax=(max(P)-min(P)),phi=phi_guess,Rd=-min(P),theta=1),trace=traceYN) 
  return(curve.nlslrc)
}

# Measured Light Intensities 
I=c(0,125,250,500,750,1000,1500,2000)
# Measured P for Plant 1
P_1=c(0.000,2.825,5.475,7.750,8.800,9.450,10.025,10.375)
# Measured P for Plant 2
P_2=c(0.050,3.000,5.100,6.050,6.350,6.675,6.850,6.925)

# Fitting LRC to Plant 1
P_Fit1 = LRcurve_predict(P_1)
print(P_Fit1$m$getPars()) #printing the parameter estimates for LRC 1

# Fitting LRC to Plant 2
P_Fit2 = LRcurve_predict(P_2,phi_guess=0.03)
print(P_Fit2$m$getPars()) # printing the parameter estimate for LRC 2


par(mar=c(3,3,0,0),oma=c(1.5,1.5,1,1))
print(plot(I,P_1,xlab="", ylab="", ylim=c(-2,15),cex.lab=1.2,cex.axis=1.5,cex=2,pch=5))
points(I,P_2,cex=2,pch=1)
lines(0:2000,predict(P_Fit1, newdata=data.frame(I=0:2000)),lwd=3,col="orange")
lines(0:2000,predict(P_Fit2, newdata=data.frame(I=0:2000)),lwd=3,col="red")
mtext(expression("Irradiance ("*mu*"mol photons "*m^-2*s^-1*")"),side=1,line=3.3,cex=2)
mtext(expression(P*" ("*mu*"mol "*CO[2]*" "*m^-2*s^-1*")"),side=2,line=2,cex=2)
title("Light Response Curves",line=-3,cex.main=2)

