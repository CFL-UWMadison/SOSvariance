# Fish Harvest, 1D alt state model, spectra
# Noise mgt case
# SRC 2014-05-14

rm(list = ls())
graphics.off()

library(multitaper)

# FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Growth function
Gfun = function(x) {
  DD = 1-(x/K) # density dependence
  con = c*x^q/(b^q + x^q) # consumption
  rate = r*x*DD - con 
  return(rate)
}

# Function to return time series given h and u
Tsim = function(h,u) {
  # Preliminaries
  eps = sigma*rnorm(nt)
  xt=rep(0,nt)
  xt[1:3]=x0 + sigma*rnorm(3)
  Ht = rep(0,nt)
  Ht[1:3] = h*xt[1:3]
  gx=rep(0,nt)
  gx[1]=Gfun(xt[1]) - Ht[1]
  gx[2]=Gfun(xt[2]) - Ht[2]
  # Simulate time series
  for( i in 3:(nt-1) )   {
    dx1=xt[i]-xt[i-1]
    dx2=xt[i-1]-xt[i-2]
    Ht[i] = h*xt[i]
    gx[i]=Gfun(xt[i]) - Ht[i]
    xnext=xt[i] + (u+v)*dx1 - u*v*dx2 + gx[i] - (u+v)*gx[i-1] + u*v*gx[i-2] + eps[i]
    xt[i+1] = max(xnext,0.1)  
  }
  # Compute last harvest
  Ht[nt] = h*xt[i-1]
  outlist = list(Ht,xt)
  return(outlist)
}

# Function to return spectra
GetSpec = function(X)  {
  Xvar = ts(X,deltat=1)
  XS = spec.mtm(Xvar, k=40, nw=20, nFFT = 'default',
                centreWithSlepians = T, Ftest = F,
                jackknife = F,maxAdaptiveIterations = 100,
                plot = F, na.action = na.fail) 
  Xspec = XS$spec 
  Xfreq=XS$freq  
  outlist = list(Xfreq,Xspec)
  return(outlist)
}

# END FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Constants for a simulation
r=0.75
c=1.35
b=1
q=2
K=8.8  # critical K is about 6.5 for h=0, 7.6 for h=0.05, 9 for h=0.1
h=0.02

w = 0 
v = 0.4
sigma = 0.1

# Start simulations

nt = (2^11) + 100
tvec=(1:(nt-100))

# Generate time series with no variance control
u = 0
x0 = (r-h)*(K/r) # deterministic equilibrium forlogistic model
simlist = Tsim(h,u)
Hsraw = simlist[[1]]
Hsim = Hsraw[101:nt] # throw away the burn-in
xsraw = simlist[[2]]
xsim = xsraw[101:nt] # throw away the burn-in
# Compute spectra for harvest
xcen = Hsim - mean(Hsim)
xstd = xcen/sd(Hsim)
# Spectral variance
sp.nocon = GetSpec(xcen)
# Save simulated harvest
Hsim.nocon=Hsim

# Generate time series with variance control
u = 0.8
x0 = (r-h)*(K/r) # deterministic equilibrium
simlist = Tsim(h,u)
Hsraw = simlist[[1]]
Hsim = Hsraw[101:nt] # throw away the burn-in
xsraw = simlist[[2]]
xsim = xsraw[101:nt] # throw away the burn-in
# Compute spectra for x
xcen = Hsim - mean(Hsim)
xstd = xcen/sd(Hsim)
# Spectral variance
sp.con = GetSpec(xcen)
# Save simulated harvest
Hsim.con=Hsim

# Generate time series with white noise
u = 0
v = 0
x0 = (r-h)*(K/r) # deterministic equilibrium
simlist = Tsim(h,u)
Hsraw = simlist[[1]]
Hsim = Hsraw[101:nt] # throw away the burn-in
xsraw = simlist[[2]]
xsim = xsraw[101:nt] # throw away the burn-in
# Compute spectra for x
xcen = Hsim - mean(Hsim)
xstd = xcen/sd(Hsim)
# Spectral variance
sp.white = GetSpec(xcen)
# Save simulated harvest
Hsim.white=Hsim

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(Hsim.nocon,Hsim.con)
plot(tvec,Hsim.nocon,type='l',lwd=2,col='blue',ylim=yrange,#log='y',
     xlab='Time Step',ylab='Harvest')
points(tvec,Hsim.con,type='l',lwd=2,col='red')
points(tvec,Hsim.white,type='l',lwd=2,col='black')


windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(sp.nocon[[2]],sp.con[[2]])
plot(sp.nocon[[1]],sp.nocon[[2]],type='l',lwd=2,col='blue',ylim=yrange,log='y',
     xlab='Frequency',ylab='Spectrum')
points(sp.con[[1]],sp.con[[2]],type='l',lwd=2,col='red')
points(sp.white[[1]],sp.white[[2]],type='l',lwd=2,col='black')
legend('topright',
       legend=c('Nominal','Var. Management','White Input'),
       col=c('blue','red','black'),
       lwd=c(2,2,2),cex=1.6)
