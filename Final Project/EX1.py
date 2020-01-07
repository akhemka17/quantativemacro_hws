
"""
Macro Final Project 

"""

import numpy as np
import quantecon
from scipy import optimize
import statsmodels.api as sm



### Paremeters
beta = (0.99)**40
alph = 0.3
lambd = 0.5
tau = 0
# tau=0.1
t=50000
mu=1
prob=np.array([0.5, 0.5]) # assuming same prob for both states

## shocks

# eta
n=11
sde=0.95
eta = quantecon.quad.qnwnorm(n, mu,sig2=sde)
probeta=eta[1]
probe=np.exp(eta[0])
# zeta
sdz=0.13
zeta = np.exp(np.random.normal(mu, sdz,2))
z=np.random.choice(zeta,t)
# rho
sdr=0.5
rho = np.exp(np.random.normal(mu, sdr,2))
r=np.random.choice(rho,t)

## Ex 1.2

# phi
phi=0
for aa, bb in enumerate(rho):
    for cc,dd in enumerate(probe):
        phi += prob[aa]*probeta[cc]*(1+(1-alph)/(alph*(1+lambd)*bb)*(lambd*dd+tau*(1+lambd*(1-dd))))**(-1)

# savings
sav = beta * phi/(1+beta*phi)

# capital in logs
logk=np.empty(t)
logk[0]=np.log(sav)+np.log(1-tau)+np.log(1-alph)
for i in range(1,t):
    logk[i] = np.log(sav)+np.log(1-tau)+np.log(1-alph)+np.log(z[i-1])+alph*(logk[i-1])


## Ex 1.3

# guess for psi
bad=np.array([zeta[0], rho[0]])
good=np.array([zeta[1], rho[1]])
z1=np.array([bad,good])
psi0=np.empty([2,2])
for i in range(0,2):
    psi0[i,0] = np.log(sav)+np.log(1-tau)+np.log(1-alph)+np.log(z1[i,0])
    psi0[i,1] = alph

## Krussel Smith Algorithm
kss=np.exp(logk[0])    
def f(psi):
    # grid space
    kgrid = np.linspace(0.5*kss,1.5*kss,5)
    k1=np.empty([2,5]) 
    savgrid=np.empty([2,5])  
    # evaluating kt+1
    for ee,ff in enumerate(kgrid):
        for jj in range(0,2):
            k1[jj,ee] = np.exp(psi[jj,0]+psi[jj,1]*np.log(ff))
    for kk in range(0,2):    
        for ee,ff in enumerate(kgrid):
            def g(ass):
                # solving for savings  
                g= 1-beta*((1-tau)*(1-alph)*ff**alph*z1[kk,0]-ass)*(1+alph*k1[kk,ee]**(alph-1)*z1[kk,0]*z1[kk,1])/(ass*(1+alph*k1[kk,ee]**(alph-1)*z1[kk,0]*z1[kk,1])+lambd*(1-alph)*k1[kk,ee]**alph*z1[kk,0]*(1-tau)+(1-lambd)*tau*(1-alph)*k1[kk,ee]**alph*z1[kk,0]*(1+lambd)/(1-lambd))
                return(g)
            savgrid[kk,ee] = optimize.brentq(g,0.0001,3)/((1-tau)*(1-alph)*ff**alph*z1[kk,0])
    return(savgrid)

# running simulation
def h(savgrid):
    sav1= np.random.choice(savgrid[0],t)
    sav2 = np.random.choice(savgrid[1],t)
    sav=np.array([sav1,sav2])
    ks=np.empty([2,t])
    c1=np.empty([2,t])
    c2=np.empty([2,t])
    ks[:,0] = kss
    for i in range(0,2):
        for j in range(1,t):
            ks[i,j] = np.exp(np.log(sav[i,j])+np.log(1-tau)+np.log(1-alph)+np.log(z[j])+alph*np.log(ks[i,j-1]))
            c1[i,j] = (1-sav[i,j])*(1-tau)*(1-alph)*z[j]*ks[i,j]**alph
            c2[i,j] = beta*c1[i,j]*(1+alph*ks[i,j]**(alph-1)*z[j]*r[j])
    return(ks[0],ks[1],c1,c2)


# regression
def i(kb,kg):
    logkb = np.log(kb)[0:t-1]
    logkg = np.log(kg)[0:t-1]
    logkg1 = np.log(kg)[1:t]
    logkb1 = np.log(kb)[1:t]
    Xg = logkg[500:]
    Xg = sm.add_constant(Xg)
    Xb = logkb[500:]
    Xb = sm.add_constant(Xb)
    Yg = logkg1[500:]
    Yb = logkb1[500:]

    mod1 = sm.OLS(Yb,Xb)
    reg1 = mod1.fit()
    mod2= sm.OLS(Yg, Xg)
    reg2 = mod2.fit()
    psig = np.array(reg2.params)
    psib = np.array(reg1.params)
    psi= np.array([psib, psig])
    return(psi)

# Solving for household    
iter= 0
maxiter=100
weight=0.7
psi1 = np.array([[3,2],[1,2]])
epsi = 0.001
psiguess=psi0

while np.max(np.abs(psi1-psiguess))>epsi:
    if iter<=maxiter:
        iter += 1
        savgrid = f(psi1)
        sol = h(savgrid)
        c1 = sol[2]
        c1b = c1[0]
        c1b = c1b[500:]
        c1g = c1[1]
        c1g = c1g[500:]
        c2 = sol[3]
        c2b = c2[0]
        c2b = c2b[500:]
        c2g = c2[1]
        c2g = c2g[500:]
        utilb = 1/(t-500)*np.sum((1/1+beta)*np.log(c1b)+(1/1+beta)*np.log(c2b))
        utilg = 1/(t-500)*np.sum((1/1+beta)*np.log(c1g)+(1/1+beta)*np.log(c2g))
        v=(utilg+utilb)/2
        psi00 = i(sol[0],sol[1])
        psiguess=psi1
        psi1 = weight*psi00+(1-weight)*psiguess
        print(iter)

# welfare
# tau=0
v0=0.960949794
v1=v
g= (v1-v0/beta)-1
        



    