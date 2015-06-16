# Program to plot spectra for 2-D lake model
# SRC 2015-05-24

rm(list = ls())
graphics.off()

library(multitaper)

# FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# 2D rate funtion for simulation 
dWM.noh = function(M,W) {
  rate = c(0,0)
  Rden=mq + W^q
  # deterministic part of dW/dt w/o load control
  rate[1] = -(s+h)*W + (r*M*W^q/Rden)
  # deterministic part of dM/dt 
  rate[2] = s*W - b*M - (r*M*(W^q)/Rden) 
  return(rate)
}

# Function to return time series given lamda and u
Tsim = function(lamda,u) {
  # Preliminaries
  eps = sigma*rnorm(nt)
  Wt=rep(0,nt)
  Mt=rep(0,nt)
  Wt[1:3]=W0 + eps[1:3]
  Mt[1:3]=M0 - eps[1:3]
  
  ratemat = matrix(0,nr=nt,ncol=2)
  ratemat[1,] = dWM.noh(Mt[1],Wt[1])
  ratemat[2,] = dWM.noh(Mt[2],Wt[2])
  
  for(i in 3:(nt-1) )  {
    ratemat[i,] = dWM.noh(Mt[i],Wt[i])
    Mnext = Mt[i] + (ratemat[i,2])*dt
    Mt[i+1] = max(Mnext,1)
    Wnext = Wt[i] + (u+phi)*(Wt[i]-Wt[i-1]) - u*phi*(Wt[i-1]-Wt[i-2]) +
      ratemat[i,1]*dt - (u + phi)*ratemat[i-1,1]*dt + u*phi*ratemat[i-2,1]*dt +
      (1-(u+phi)+u*phi)*lamda + eps[i]
    Wt[i+1] = max(Wnext,0.1)
  }
  outlist=list(Wt,Mt)
  return(outlist)
}

# Function to return spectra
GetSpec = function(X)  {
  Xvar = ts(X,deltat=1)
  XS = spec.mtm(Xvar, k=30, nw=15, nFFT = 'default',
                centreWithSlepians = T, Ftest = F,
                jackknife = F,maxAdaptiveIterations = 100,
                plot = F, na.action = na.fail) 
  Xspec = XS$spec 
  Xfreq=XS$freq  
  outlist = list(Xfreq,Xspec)
  return(outlist)
}

# END FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Parameters 
b = 0.002 # this rate is increased from previous papers; burial is 0.001 in E.L. 2006
#h = 0.15 # export coef. from E.L. 2006
h = 1-exp(-0.29) # export coef. from BalanceModelFit_Mendota+2sigmas_2013-08-10.R (Stable balance model)
m = 4 # half-saturation for recycle; wide range 1-15 in Carpenter & Lathrop Ecosystems 2008; 8 in EL 2006
r = 0.019 # max. recycling rate; ~0.002 in Carpenter & Lathrop 2008; 0.019 in E.L. 2006
q = 4 # 4 is near posterior mode in Carpenter & Lathrop Ecosystems 2008; 8 was used in E.L. 2006
mq = m^q  
#s = 0.7  # sedimentation from E.L. 2006
#s = 1-exp(-0.34) # sedimentation from BalanceModelFit_Mendota+2sigmas_2013-08-10.R (Stable balance model)
s = 1-h # see typed notes on this model from 2014-12-06

# Variables chosen to be near the threshold but within safe operating space
# See Lake2D_U-in-h_ARMA_V0_2014-12-23.R
# Note that mean Mendota load is about 0.85 g m-2 y-1 (range 0.32 to 1.98)
L.mod = 1.2
M.mod = 330

# Noise process
sigma = 0.35
phi = 0.1 # AR

# Simulation control
nt = (2^11)+500
dt=1
tvec = (1:(nt-500))

# Initial conditions
W0 = 1
M0 = M.mod

# Set up lamda value
lamda=0.5

# Set up u values
nu=2
uvec = c(0.6,0)

# Compute output statistics 
Wmat = matrix(0,nr=nt-500,nc=2)
Mmat = matrix(0,nr=nt-500,nc=2)

# Compute time series
  for(j in 1:nu)  { # loop over u
  u=uvec[j]
  simlist = Tsim(lamda,u)
  Wsraw = simlist[[1]]
  Wmat[,j] = Wsraw[501:nt] # throw away the burn-in
  Msraw = simlist[[2]]
  Mmat[,j] = Msraw[501:nt] # throw away the burn-in
  } # end loop over u

# Compute spectra
# Variance control
xcen = Wmat[,1] - mean(Wmat[,1])
xstd = xcen/sd(Wmat[,1])
sp.con = GetSpec(xstd)
# No Variance control
xcen = Wmat[,2] - mean(Wmat[,2])
xstd = xcen/sd(Wmat[,2])
sp.nocon = GetSpec(xstd)

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(Wmat)
plot(tvec,Wmat[,2],type='l',lwd=2,col='blue',ylim=yrange,#log='y',
     xlab='Time Step',ylab='Water P')
points(tvec,Wmat[,1],type='l',lwd=2,col='red')
#points(tvec,Hsim.white,type='l',lwd=2,col='black')

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(sp.nocon[[2]],sp.con[[2]])
plot(sp.nocon[[1]],sp.nocon[[2]],type='l',lwd=2,col='blue',ylim=yrange,log='y',
     xlab='Frequency',ylab='Spectrum')
points(sp.con[[1]],sp.con[[2]],type='l',lwd=2,col='red')
legend('bottomleft',
       legend=c('Nominal','Var. Management'),
       col=c('blue','red'),
       lwd=c(2,2),cex=1.6)
