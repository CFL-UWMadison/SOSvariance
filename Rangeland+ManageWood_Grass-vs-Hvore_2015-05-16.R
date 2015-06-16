# Program to plot crashplots for Rangeland model
# This version manages the wood, e.g. by herbicide or cutting
# SRC 2015-05-16

rm(list = ls())
graphics.off()

library(multitaper)

# FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Grass and wood dynamics with herbivory and wood mortality omitted
dGW.noH = function(G,W) {
  # infiltration
  I = (G + c1*i0)/(G + c1)
  # grass
  UG = 1/(G + alpha*W + beta)
  Qnz = c2*(1 - (G0/G))/(G + c3)
  Q = max(Qnz,0)
  # wood
  S = rS*( 1 - (G/(K*(G + d4))) )
  UW = 1/(W  + gamma*G + delta)
  # rates
  dGdt = rG*I*UG*G - LG*G + Gin #- H*Q*G 
  dWdt = rW*I*UW*W + (S*W/(W+p)) #- LW*W
  rates = c(dGdt,dWdt,Q)
  return(rates)
}

# Function to return time series given h and u
Tsim = function(H,u) {
  uG = 0 # grass u is zero for this version of the model
  # Preliminaries
  eps = sigma*rnorm(nt)
  Gt=rep(0,nt)
  Wt=rep(0,nt)
  Hvory = rep(0,nt)
  Qt = rep(0,nt)
  Gt[1:3]=Gt0 + eps[1:3]
  Wt[1:3]=Wt0 
  ratemat = matrix(0,nr=nt,ncol=2)
  rates = dGW.noH(Gt[1],Wt[1])
  ratemat[1,] = rates[1:2]
  Qt[1] = rates[3]
  Hvory[1] = H*Gt[1]*Qt[1]
  rates = dGW.noH(Gt[2],Wt[2])
  ratemat[2,] = rates[1:2]
  Qt[2] = rates[3]
  Hvory[2] = H*Gt[2]*Qt[2]
  for(i in 3:(nt-1) )  {
    rates = dGW.noH(Gt[i],Wt[i])
    ratemat[i,] = rates[1:2]
    Qt[i] = rates[3]
    Hflux = H*Qt[i]*Gt[i] - (uG+phi)*H*Qt[i-1]*Gt[i-1] + uG*phi*H*Qt[i-2]*Gt[i-2]
    Hvory[i+1] = max(Hflux,0.1)
    Wnext = Wt[i] + (ratemat[i,2])*dt - LW*Wt[i]*dt + u*LW*Wt[i-1]
    Wt[i+1] = max(Wnext,1)
    Gnext = Gt[i] + phi*(Gt[i]-Gt[i-1]) + (ratemat[i,1]*dt) - phi*ratemat[(i-1),1]*dt - 
      Hvory[i]*dt + eps[i] - theta*eps[i-1]
    Gt[i+1] = max(Gnext,1)
  }
  outlist=list(Gt,Wt,Hvory)
  return(outlist)
}

# Function to return spectra
GetSpec = function(X)  {
  Xvar = ts(X,deltat=1)
  XS = spec.mtm(Xvar, k=15, nw=30, nFFT = 'default',
                centreWithSlepians = T, Ftest = F,
                jackknife = F,maxAdaptiveIterations = 100,
                plot = F, na.action = na.fail) 
  Xspec = XS$spec 
  Xfreq=XS$freq  
  outlist = list(Xfreq,Xspec)
  return(outlist)
}

# END FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Main program *****************************************************

# Parameters from Walker et al. 1981
# Most are found in Table 1

# grass
rG = 4000 # growth coef. from topsoil water
c2 = 3.15 # consumption rate by herbivore
G0 = 500  # grazing threshold
d3 = 600 # half-saturation G for consumption
c3 = d3 - G0 # p. 485 under eq 3
i0 = 0.1 # minimum infiltration at low values of G
alpha = 0.122/0.21  # water uptake efficiency by wood / w.u.e. of grass
LG = 1 # annual grass loss by respiration and death
# Beta is e / theta_G where e is rate of evaporation of soil
# water (as depth/time) and theta_G is 0.21 according to Table 1
# On page 483 under eq. 2 they state "evaporation loss is generally
# insignificant . . . "
beta = 0.01/0.21  
# c1, the rate at which infiltration increases with G
d1 = 500
d2 = i0*d1/(1-i0)
c1 = d1+d2
# Grass input (colonization) NOT in original model
Gin = 1

# wood
rW = 3000  # growth coef. from topsoil water
rS = 2000  # growth coef. from subsoil water (4000 in Walker et al)
d4 = 500  # half saturation effect of G on S
K = 2  # proportional reduction of rS at high G
p = 0.1 # inefficiency of subsoil water use
LW = 1  # annual loss of wood veg by resp & death
gamma = 1/alpha # top line p. 487
delta = beta/alpha # top line p. 487

# Noise process
sigma = 5
theta = 0 # MA
phi = 0.2 # AR

# Simulation control
nt=1500
dt=1

# Initial conditions
Gt0 = 2000
Wt0 = 2000

# Set up H values
nH = 15
Hvec = seq(1,100,length.out=nH)

# Set up u values
nu=3
uvec = c(0.07,0,-0.5)

# Compute output statistics 
Hh = matrix(0,nr=nH,nc=nu) # mean Herbivore flux 
Tvar = matrix(0,nr=nH,nc=nu)  # total variance of Herbivore flux
HFvar = matrix(0,nr=nH,nc=nu)  # high freq variance of herbivore flux
HFutil = matrix(0,nr=nH,nc=nu) # utility based on high freq var at each h
Gmean = matrix(0,nr=nH,nc=nu)
Wmean = matrix(0,nr=nH,nc=nu)
PhiG = matrix(0,nr=nH,nc=nu) # proportion of days with grass > G0

for(i in 1:nH) { # loop over h
  for(j in 1:nu)  { # loop over u
  H=Hvec[i]
  u=uvec[j]
  simlist = Tsim(H,u)
  Gsraw = simlist[[1]]
  Gsim = Gsraw[501:nt] # throw away the burn-in
  Wsraw = simlist[[2]]
  Wsim = Wsraw[501:nt] # throw away the burn-in
  Hsraw = simlist[[3]]
  Hsim = Hsraw[501:nt] # throw away the burn-in
  Hh[i,j] = mean(Hsim)
  Gmean[i,j] = mean(Gsim)
  Wmean[i,j] = mean(Wsim)
  NGgood = ifelse(Gsim>G0,1,0)
  PhiG[i,j] = sum(NGgood)/(nt-500)
  # Compute spectra
  xcen = Hsim - mean(Hsim)
  xstd = xcen/sd(Hsim)
  # Spectral variance
  sp = GetSpec(xcen)
  Tvar[i,j] = sum(sp[[2]])
  sp.HF = subset(sp[[2]],subset=(sp[[1]]<=0.5 & sp[[1]]>=0.25))
  HFvar[i,j] = sum(sp.HF)
  }  # end loop over u
} # end loop over h

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(Gmean)
plot(Hvec,Gmean[,1],type='b',lwd=2,pch=19,col='red',ylim=yrange,
     xlab='Stocking Level',ylab='Grass Biomass')
points(Hvec,Gmean[,2],type='b',lwd=2,pch=21,cex=1.3,col='blue')
points(Hvec,Gmean[,3],type='b',lwd=2,pch=17,cex=1.3,col='black')

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(Wmean)
plot(Hvec,Wmean[,1],type='b',lwd=2,pch=19,col='red',ylim=yrange,#log='x',
     xlab='Stocking Level',ylab='Woody Vegetation')
points(Hvec,Wmean[,2],type='b',lwd=2,pch=19,col='blue')