#!/usr/bin/env python2

import numpy as np
import sys
import scipy.constants
import os

def mean(a,weights):
  return np.average(a,weights=weights)

def meansq(a,weights):
  return np.average(a**2,weights=weights)

def var(a,weights):
  return np.average(a**2,weights=weights)-np.average(a,weights=weights)**2

def rms(a,weights):
  return np.sqrt(meansq(a,weights=weights))

def std(a,weights):
  return np.sqrt(var(a,weights=weights))


phase=float(sys.argv[1])
if phase<0:
  filename="beamfile.bin"
else:
  filename="tb%s.swp" % sys.argv[1].zfill(5)


c=scipy.constants.c
e=scipy.constants.e
m=scipy.constants.electron_mass
mp=scipy.constants.proton_mass
e0=scipy.constants.epsilon_0
pi=np.pi
IA=4*pi*e0*m*c**3/e

if 'LCODE_ne' in os.environ:
  n=float(os.environ['LCODE_ne'])
else:
  n=2e20 #/m^3
if 'LCODE_dz' in os.environ:
  dz=float(os.environ['LCODE_dz'])
else:
  dz=0.02

op=(n*e**2/e0/m)**.5
k=op/c
um=1e6/k

z,r,pz,pr,pa,q,w,N=np.fromfile(filename).reshape(-1,8)[:-1].transpose()
#print (np.mean(r**2,weight=mw)*np.var(pr)-np.mean(r*(pr-np.mean(pr)))**2)**.5/1836/2**.5*200.88

zbins=np.linspace(-80,0,400)
Nbins=len(zbins)

nbin=np.digitize(z,zbins)

print(max(z), min(z))
print(max(zbins), min(zbins))

binned=np.zeros((Nbins,5)) # rms(r), em, alpha, w, N

print(nbin)

#np.savetxt("test.dat",np.c_[nbin])

#sys.exit(0)

sq=-1
for bin in range(Nbins):
  print("slice",bin,"of",Nbins)
  mw=w*(q==sq)*(nbin==bin)*IA*dz/2/op*1e12 # masked weight
  N=np.count_nonzero(mw)
  
  if N==0:
    continue
  
  msr =meansq(r,weights=mw)*um**2
  mspr=meansq(pr,weights=mw)
  mspa=meansq(pa/r,weights=mw)

  
  am = -mean(r*pr,weights=mw)
  em  =(msr*(mspr+mspa)-am**2)**.5/2 # from LCODE manual
  alpha = -mean(r*pr,weights=mw)/em
  
  if(N):
    binned[bin,:]=[msr**.5/2**.5,\
                  em, alpha,
                  np.sum(mw),\
                  N]
  
np.savetxt("ems%d.dat" % phase,np.c_[zbins,binned],delimiter=" ")

print("Weights: %f %f" % (np.sum(w)*IA*dz/2/op*1e12, np.sum(binned[:,3])))
