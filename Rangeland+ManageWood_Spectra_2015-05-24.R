# Program to plot spectra for Rangeland model
# Woody is managed, e.g. by herbicide or cutting
# SRC 2015-05-24

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
nt = (2^11)+500
dt=1
tvec = (1:(2^11))

# Initial conditions
Gt0 = 2000
Wt0 = 2000

# Set up H value
H = 10

# Set up u values
nu=3
uvec = c(0.07,0,-0.5)

# Compute output statistics 
Hmat = matrix(0,nr=2^11,nc=nu)
Gmat = matrix(0,nr=2^11,nc=nu)
Wmat = matrix(0,nr=2^11,nc=nu)

# Compute time series for analysis
  for(j in 1:nu)  { # loop over u
  u=uvec[j]
  simlist = Tsim(H,u)
  Gsraw = simlist[[1]]
  Gsim = Gsraw[501:nt] # throw away the burn-in
  Wsraw = simlist[[2]]
  Wsim = Wsraw[501:nt] # throw away the burn-in
  Hsraw = simlist[[3]]
  Hsim = Hsraw[501:nt] # throw away the burn-in
  Hmat[,j] = Hsim
  Gmat[,j] = Gsim
  Wmat[,j] = Wsim
  }  # end loop over u


# Compute spectra
# Choose output matrix
Xmat = Gmat
# Variance control
xcen = Xmat[,1] - mean(Xmat[,1])
xstd = xcen/sd(Xmat[,1])
sp.con = GetSpec(xstd)
# No Variance control
xcen = Xmat[,2] - mean(Xmat[,2])
xstd = xcen/sd(Xmat[,2])
sp.nocon = GetSpec(xstd)
# Whitened spectrum
xcen = Xmat[,3] - mean(Xmat[,3])
xstd = xcen/sd(Xmat[,3])
sp.white = GetSpec(xstd)

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(Xmat)
plot(tvec,Xmat[,2],type='l',lwd=2,col='blue',ylim=yrange,#log='y',
     xlab='Time Step',ylab='Output Variable')
points(tvec,Xmat[,1],type='l',lwd=2,col='red')
points(tvec,Xmat[,3],type='l',lwd=2,col='black')

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange=range(sp.nocon[[2]],sp.con[[2]])
plot(sp.nocon[[1]],sp.nocon[[2]],type='l',lwd=2,col='blue',ylim=yrange,log='y',
     xlab='Frequency',ylab='Spectrum')
points(sp.con[[1]],sp.con[[2]],type='l',lwd=2,col='red')
points(sp.white[[1]],sp.white[[2]],type='l',lwd=2,col='black')
legend('topright',
       legend=c('Nominal','u = 0.07','u = -0.5'),
       col=c('blue','red','black'),
       lwd=c(2,2),cex=1.6)
