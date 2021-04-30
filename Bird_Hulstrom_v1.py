import numpy as np
import datetime as dat


n=dat.date(2021,2,9).timetuple().tm_yday 
h=0
Wo=0.9
T=305.15
HR=0.59
alfa=1.3
beta=0.4
zen=(np.pi/2) - np.deg2rad(65.14)
Lo=239/100000 #This is either /100 or /100000
rhog=0.02


Cr=1367*(1+0.033*np.cos(n*360/365))

mrel=1/(np.cos(zen)+0.15*(93.885-zen)**(-1.253))
Pt=101325*np.exp(-0.0001184*h)
ma=Pt*mrel/101325
WW=0.493*(HR/T)*np.exp(26.23-5416/T)
Uw=WW*mrel
Fc=0.93-0.21*np.log(ma)

Tr=np.exp(-0.0903*(ma**0.84)*(1+ma-ma**1.01))
To=1 - 0.1611*((Lo*mrel)*(1+139.48*Lo*mrel)**(-0.3035)) - (0.002715*Lo*mrel)/(1+0.044*Lo*mrel+0.003*(Lo*mrel)**2)
Tg=np.exp(-0.0127*ma**0.26)
Tw=1- 2.4959*Uw/((1+79.034*Uw)**0.6828 + 6.385*Uw)
Ta=0.12445*alfa-0.0162 + (1.003-0.125*alfa)*np.exp(-beta*ma*(1.089*alfa+0.5123))

Idh=0.9662*Cr*Tr*To*Tg*Tw*Ta*np.cos(zen)



Taa=1-(1-Wo)*(1-ma+ma**1.06)*(1-Ta)
Tas=Ta/Taa
rho_alfa_p=0.0685+(1-Fc)*(1-Tas)


Idr=0.79*Cr*To*Tg*Tw*Taa*0.5*((1-Tr)/(1-ma+ma**1.02))*np.cos(zen)
Ida=0.79*Cr*To*Tg*Tw*Taa*Fc*((1-Tas)/(1-ma+ma**1.02))*np.cos(zen)
Idm=(Idh*np.cos(zen) *Idr+ Ida)*(rhog*rho_alfa_p/(1-rhog*rho_alfa_p))

Idfh=Idr+Ida+Idm