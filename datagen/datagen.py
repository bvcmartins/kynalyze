import sys
import os
import numpy as np
import random
import pylab as pb
#import scipy.fftpack
#import scipy.special
#from scipy.optimize import curve_fit
#from scipy.optimize import leastsq
#import pickle
#import matplotlib.pyplot as plt

'''
This should generate a bunch of artificial data, and write it to the folder ./data. This 
does not assume any particular number of states.
'''

def sigofmu(mu):
  sig=0.5+(1.0/20)*mu
  return sig

def gendat(P,mu,sig,GAMMAS,filename):
  Npts=20000
  Fs=10000
  Nstates=len(mu)
  test=random.random()
  thresh=[];Pcumul=0.0
  # set initial state
  for k in range(Nstates): thresh.append(Pcumul) ; Pcumul+=P[k]
  for k in range(Nstates):
    if test>thresh[k]: state=k
  
  I=[]
  ST=[]
  for k in range(Npts):
    y=random.gauss(mu[state],sig[state])
    ST.append(mu[state])
    I.append(y)
    for j in range(Nstates):
      if random.random()<GAMMAS[state][j]/Fs: state=j
  f=open(filename,'w')
  tau=0.0
  for i in I:
    f.write(str(tau)+' '+str(i*1e-12)+'\n')
    tau+=1.0/Fs
  f.close()
  #pb.plot(I,'b')
  #pb.plot(ST,'r')
  #pb.xlim([0,1000])
  #pb.show()
  return
  

def main():
  MINDIST=2.0 #nm
  MAXDIST=6.0 #nm
  NUMPTS=60
  basename='./fakedata/Agen'
  f=open(basename+'.dat','w')
  f.write("Starting Point:"+'      '+str(MINDIST)+'      0.0\n')
  f.write("Ending Point:"+'      '+str(MAXDIST)+'      0.0\n')
  f.write("Number of Points:"+'   '+str(NUMPTS)+'\n')
  f.write("Voltage:"+' '+'0.0')
  f.close()
  q=0;
  for DIST in np.linspace(MINDIST,MAXDIST,NUMPTS):
    filename=basename+"%03d"%(q+1)+'.dat'
    print ''
    print filename
    
    #define parameters
    mu1=11.0+DIST*10.0/8.0
    mu2=21.0
    mu3=41.0-DIST*20.0/8.0
    mu=(mu1,mu2,mu3) # no particular reason for these values
    sig=tuple([sigofmu(x) for x in mu])
    GAMMAS=np.zeros((len(mu),len(mu))) #the values will be the rates
    FILL=500000.0*np.exp(-DIST/0.5)
    EMPTY=500.0
    for i in range(len(mu)):
      for j in range(len(mu)):
        if i==j: GAMMAS[i][j]=0.0
        if (i-j)==1: GAMMAS[i][j]=FILL # transitions between adjacent states
        if (i-j)==-1: GAMMAS[i][j]=EMPTY # transitions between adjacent states
    #we should get P by getting M from GAMMAS then doing linear algebra, but for now...
    M=np.copy(GAMMAS)
    for i in range(len(mu)):
      TEMP=0.0
      for k in range(len(mu)):
        TEMP+=GAMMAS[k][i]
      M[i][i]=-TEMP
    J,S=np.linalg.eig(M)
    J=list(np.abs(J))
    Pindex=J.index(min(J))
    P=[]
    for k in range(len(mu)):
      P.append(S[k][Pindex])
    P.reverse() # this seems to make things work but I don't know why.
    P=np.array(P)/sum(P)
    #P=tuple(np.ones(len(mu),dtype=float)/len(mu)) # equal probabilities of each state
  
    #check that parameters are reasonable
    if np.sum(P)<0.99 or np.sum(P)>1.01: print("probabilities don't add to 1"); sys.exit(0)
    if len(mu)!=len(P) or len(P)!=len(sig): print("P and mu and sig need to be the same length"); sys.exit(0)
  
    #generate artificial data
    gendat(P,mu,sig,GAMMAS,filename)
    f=open(basename+"%03d"%(q+1)+'_ANSWER.dat','w')
    f.write("FILL: "+str(FILL)+'\n')
    f.write("EMPTY: "+str(EMPTY)+'\n')
    f.write("the means are:\n")
    for mui in mu:
      f.write("  "+str(mui)+'\n')
    f.close()
    q+=1
  return

if __name__=='__main__':
  main()
