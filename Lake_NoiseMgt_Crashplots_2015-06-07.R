# Program for crash plots for 2-D lake model
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
  XS = spec.mtm(Xvar, k=7, nw=4.0, nFFT = 'default',
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
M.mod = 300

# Noise process
sigma = 0.35
phi = 0.1 # AR

# Simulation control
nt=5500
dt=1
tvec = (1:(nt-500))

# Initial conditions
W0 = 1
M0 = M.mod

# Set up lamda values
nlam = 15
lamvec = seq(0.25,2.5,length.out=nlam)

# Set up u values
# Broad search
nu=2
uvec = c(0.6,0)

# Compute output statistics 
arisk = 0.5 # risk aversion coefficient
Hh = matrix(0,nr=nlam,nc=nu) # mean water P at each h
Tvar = matrix(0,nr=nlam,nc=nu)  # total variance at each h
HFvar = matrix(0,nr=nlam,nc=nu)  # high freq variance at each h
Pctoligo = matrix(0,nr=nlam,nc=nu) # proportion of oligotrophic days

for(i in 1:nlam) { # loop over lamda
  for(j in 1:nu)  { # loop over u
  lamda = lamvec[i]
  u=uvec[j]
  simlist = Tsim(lamda,u)
  Wsraw = simlist[[1]]
  Wsim = Wsraw[501:nt] # throw away the burn-in
  Msraw = simlist[[2]]
  Msim = Msraw[501:nt] # throw away the burn-in
  Hh[i,j] = mean(Wsim)
  Noligo = ifelse(Wsim<=3,1,0)
  Pctoligo[i,j] = sum(Noligo)/(nt-500)
  # Compute spectra
  xcen = Wsim - mean(Wsim)
  xstd = xcen/sd(Wsim)
  # Spectral variance
  sp = GetSpec(xstd)
  Tvar[i,j] = sum(sp[[2]])
  sp.HF = subset(sp[[2]],subset=(sp[[1]]<=0.5 & sp[[1]]>=0.25))
  HFvar[i,j] = sum(sp.HF)
  }  # end loop over u
} # end loop over lamda

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(lamvec,Pctoligo[,1],type='b',lwd=2,pch=19,col='red',
     xlab='P Load Rate',ylab='Proportion Oligotrophic')
points(lamvec,Pctoligo[,2],type='b',lwd=2,pch=21,cex=1.3,col='blue')
#legend('bottomleft',
#       legend=c('Nominal','Var. Management'),
#       col=c('blue','red'),
#       lwd=c(2,2),cex=1)
