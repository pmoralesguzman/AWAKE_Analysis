#!/usr/bin/env python

"""
   Generate the Run2 toy model beam distributions for lcode

   Written by John.
"""

import numpy as np
import scipy
import scipy.constants
import scipy.stats as stats
import sys

c=scipy.constants.c
e=scipy.constants.e
m=scipy.constants.electron_mass
mp=scipy.constants.proton_mass
e0=scipy.constants.epsilon_0
pi=np.pi
IA=4*pi*e0*m*c**3/e
endparticle=[-100000.0,0,0,0,0,1,0,0]

n=7e20 #/m^3
op=(n*e**2/e0/m)**.5
k=op/c

dz=0.02
dr=0.02

if len(sys.argv) < 3:
  print "not enough inputs"
  print sys.exit(0)

Lw=60e-6
#Rw=10e-6*2**.5
z0w=-1 #1.25-Lwindow+Lbase
Ew=150e6*e/m/c**2
emw=float(sys.argv[1])*1e-6    #emittance in um  #6.84e-6
Qw =float(sys.argv[2])*-1e-12  #charge in pC     #-100e-12
Fpos=float(sys.argv[3])

Nw=100000
#Rw=8.66969*emw**.5/Ew**.25/(n/1e6)**.25
Rw=(2*c**2*emw**2/op**2/Ew)**.25
print "Matched radius :",Rw*1e6,"um"

zw=stats.norm.rvs(scale=Lw*k,loc=z0w,size=Nw)
zw=np.sort(zw)[::-1]
#rw=stats.weibull_min.rvs(2,scale=2**.5*Rw*k,size=Nw)
xw=stats.norm.rvs(scale=Rw*k,size=Nw)
yw=stats.norm.rvs(scale=Rw*k,size=Nw)

pzw=Ew*np.random.normal(loc=1,scale=0.001,size=Nw)
#prw=emw/Rw*np.random.normal(size=Nw)
#paw=emw/Rw*np.random.normal(size=Nw)*rw
pxw=emw/Rw*np.random.normal(size=Nw)
pyw=emw/Rw*np.random.normal(size=Nw)

xw+=pxw/pzw*Fpos
yw+=pyw/pzw*Fpos
rw=(xw**2+yw**2)**.5
prw=(xw*pxw+yw*pyw)/rw
paw=(xw*pyw-yw*pxw)

qw=-np.ones(Nw)
print "Normalised particle weight =",2*Qw*c/IA / (dz/k) / Nw,"x",Nw
ww=2*Qw*c/IA / (dz/k) / Nw * np.ones(Nw)


witness=np.c_[zw,rw,pzw,prw,paw,qw,ww,range(1,Nw+1)]
endparticle=[-100000.0,0,0,0,0,1,0,0]
particles=np.r_[witness,[endparticle]]

particles.tofile("beamfile.bin")
with open("beamfile.bit",'w') as f:
  f.write("0.0")
