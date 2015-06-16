# Program to measure time to collapse for Fish Harvest 1D alt state model
# Harvest mgt case
# SRC 2015-05-14

rm(list = ls())
graphics.off()

# FUNCTIONS <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Growth function
Gfun = function(x) {
  DD = 1-(x/K) # density dependence
  con = c*x^q/(b^q + x^q) # consumption
  rate = r*x*DD - con 
  return(rate)
}

# Function to count the number of crashes in a time series
CrashCount = function(h,u) {
  Cvec = rep(0,nt)
  # Preliminaries
  eps = sigma*rnorm(nt)
  xt=rep(0,nt)
  xt[1:3]=x0 + sigma*rnorm(3)
  Ht = rep(0,nt)
  Ht[1:3] = h*xt[1:3]
  gx=rep(0,nt)
  gx[1]=Gfun(xt[1])
  gx[2]=Gfun(xt[2])
  # Simulate time series
  for( i in 3:(nt-1) )   {
    dx1=xt[i]-xt[i-1]
    dx2=xt[i-1]-xt[i-2]
    gx[i]=Gfun(xt[i])
    Ht[i+1] = -1*(- h*xt[i] + (u+v)*h*xt[i-1] - u*v*h*xt[i-2])
    xnext=xt[i] + v*dx1 + gx[i] - v*gx[i-1] - h*xt[i] + (u+v)*h*xt[i-1] - u*v*h*xt[i-2] + eps[i] - w*eps[i-1] 
    Cvec[i] = ifelse(xnext<=0.12,1,0)
    xt[i+1] = ifelse(xnext<=0.12,x0,xnext)  
  }
  Ncrash = sum(Cvec)
  return(Ncrash)
}

# END FUNCTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Constants for a simulation
r=0.75
c=1.35
b=1
q=2
K=8.8  # critical K is about 6.5 for h=0, 7.6 for h=0.05, 9 for h=0.1
h=0.06

x0 = (r-h)*(K/r) # deterministic equilibrium for logistic model

w = 0 
v = 0.4
u = 0
sigma = 0.1

# Count crashes over a vector of h
nh = 10
hvec = seq(0.01,0.1,length.out=nh)
nt = 100000 

# Baseline case, no control
Crash.nocon = rep(0,nh)

for(i in 1:nh) {
  h = hvec[i]
  Crash.nocon[i] = CrashCount(h,u)
}

T2crash.nocon = nt/(Crash.nocon+1)

# Variance-controlled case
u = -0.8
Crash.con = rep(0,nh)

for(i in 1:nh) {
  h = hvec[i]
  Crash.con[i] = CrashCount(h,u)
}

T2crash.con = nt/(Crash.con+1)

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
yrange = range(0,T2crash.nocon,T2crash.con)
plot(hvec,T2crash.con,type='b',pch=19,cex=1,lwd=2,col='red',ylim=yrange,xlim=c(0,0.1),
     xlab='h, harvest coefficient',
     ylab='Average Time Steps to Collapse')
points(hvec,T2crash.nocon,type='b',pch=19,cex=1,lwd=2,col='blue')


