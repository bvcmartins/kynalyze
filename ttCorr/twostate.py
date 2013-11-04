import sys
import os
import numpy as np
import pylab as pb
import scipy.fftpack
import scipy.special
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
#import pickle
import matplotlib.pyplot as plt

from ttCorr import *

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

numbin=300
Imax=100
Fs=1e4

def adjustGuess(x,n,guessparams):
  print guessparams
  good=False
  while good==False:
    var=raw_input("Are you pleased with the guess ('y', 'show',  or 'parameter=value') ? ")
    if 'A1=' in var: 
      guessparams=list(guessparams)
      guessparams[0]=float(var[3:])
      guessparams=tuple(guessparams)
    elif 'A2=' in var: 
      guessparams=list(guessparams)
      guessparams[2]=float(var[3:])
      guessparams=tuple(guessparams)
    elif 'mu1=' in var: 
      guessparams=list(guessparams)
      guessparams[1]=float(var[4:])
      guessparams=tuple(guessparams)
    elif 'mu2=' in var: 
      guessparams=list(guessparams)
      guessparams[3]=float(var[4:])
      guessparams=tuple(guessparams)
    elif var=='show':
      showHist(x,n,tuple(guessparams)+(0,-1000))
    elif var=='y':
      good=True
    elif var!='':
      print "Use one of the accepted parameters (A1,mu1,A2,mu2)"
  return guessparams

#def showIrange(x,n,params,I,Iranges):
#  IrangeNEG,IrangeNEU,IrangePOS=Iranges
#  pb.clf()
#  pb.bar(x,n,float(Imax)/numbin,color='0.8',linewidth=0.4,align='center')
#  nCorrNEG,binsCorr=corrHist(I,IrangeNEG,0)
#  nCorrNEU,binsCorr=corrHist(I,IrangeNEU,0)
#  nCorrPOS,binsCorr=corrHist(I,IrangePOS,0)
#  pb.bar(x,nCorrNEG,float(Imax)/numbin,color=col1,linewidth=0.4,align='center')
#  pb.bar(x,nCorrNEU,float(Imax)/numbin,color=col2,linewidth=0.4,align='center')
#  pb.bar(x,nCorrPOS,float(Imax)/numbin,color=col3,linewidth=0.4,align='center')
#  pb.plot(np.linspace(0,40,500),oneGauss(params[0:3],np.linspace(0,40,500)),color=col1,linewidth=1)
#  pb.plot(np.linspace(0,40,500),oneGauss(params[3:6],np.linspace(0,40,500)),color=col2,linewidth=1)
#  pb.plot(np.linspace(0,40,500),oneGauss(params[6:9],np.linspace(0,40,500)),color=col3,linewidth=1)
#  pb.plot(np.linspace(0,40,500),threeGauss(params,np.linspace(0,40,500)),'--k',linewidth=1)
#  pb.xlim([0,40])
#  pb.grid(True)
#  pb.show()
#  return

def getCorr(I,Irange,taumax,params):
  q=(params[1],params[3])
  n,bins=getHist(I)
  x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
  params=twoGaussFit(x,n)
  T,P1,P2=([],[],[])
  for tau in range(taumax+1)[4:]:
    print '\r'+str('%d%%'%(100.*tau/taumax)),
    sys.stdout.flush()
    nCorr,binsCorr=corrHist(I,Irange,tau)
    params2=twoGaussFit_fixed(x,nCorr,q)
    T.append(tau)
    P1.append(params2[0]/sum(params2))
    P2.append(params2[1]/sum(params2))
  return (T,P1,P2)
  
# NOT USED FOR NOW (make a mode for this) 
#def threePLOTS_corrHist(dirname,I,Irange1,Irange2,Irange3,taumax):
#  if dirname[-1]!='/':
#    dirname=dirname+'/'
#  for tau in range(taumax+1):
#    print 'working on tau='+str('%d'%tau)
#    fname=dirname+str('%03d'%tau)+'.png'
#    n,bins=getHist(I)
#    pb.clf()
#    pb.bar(bins[1:],n,float(Imax)/numbin,color='0.5',alpha=0.5,linewidth=0)
#    nCorr1,binsCorr1=corrHist(I,Irange1,tau)
#    nCorr2,binsCorr2=corrHist(I,Irange2,tau)
#    nCorr3,binsCorr3=corrHist(I,Irange3,tau)
#    pb.bar(binsCorr1[1:],nCorr1,float(Imax)/numbin,color=col1,linewidth=0.4,alpha=0.5)
#    pb.bar(binsCorr2[1:],nCorr2,float(Imax)/numbin,color=col2,linewidth=0.4,alpha=0.5)
#    pb.bar(binsCorr3[1:],nCorr3,float(Imax)/numbin,color=col3,linewidth=0.4,alpha=0.5)
#    pb.grid(True)
#    pb.xlabel('I (pA)',fontsize=16)
#    pb.xlabel('counts',fontsize=16)
#    plt.text(33,0.9*max(n),'tau='+str('%d'%tau))
#    pb.savefig(fname)
#  return

def ranges(params):
  A1,mu1,A2,mu2=params
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  IrangeONE=(0,mu2-2*sig2)
  IrangeTWO=(mu2-sig2,mu2+6*sig2)
  #ONE integrals
  ONEone=A1*(scipy.special.erf((IrangeONE[1]-mu1)/(np.sqrt(2)*sig1))-scipy.special.erf((IrangeONE[0]-mu1)/(np.sqrt(2)*sig1)))
  ONEtwo=A2*(scipy.special.erf((IrangeONE[1]-mu2)/(np.sqrt(2)*sig2))-scipy.special.erf((IrangeONE[0]-mu2)/(np.sqrt(2)*sig2)))
  tot=ONEone+ONEtwo
  PoONE=[ONEone/tot,ONEtwo/tot]#can add a third entry of zero perhaps
  #TWO integrals
  TWOone=A1*(scipy.special.erf((IrangeTWO[1]-mu1)/(np.sqrt(2)*sig1))-scipy.special.erf((IrangeTWO[0]-mu1)/(np.sqrt(2)*sig1)))
  TWOtwo=A2*(scipy.special.erf((IrangeTWO[1]-mu2)/(np.sqrt(2)*sig2))-scipy.special.erf((IrangeTWO[0]-mu2)/(np.sqrt(2)*sig2)))
  tot=TWOone+TWOtwo
  PoTWO=[TWOone/tot,TWOtwo/tot]
  
  PoALL=(PoONE,PoTWO)
  return (IrangeONE,IrangeTWO,PoALL)

def singlexPLOT(ONE,TWO,p,q,PoALL,taumax,directory=None):
  #ONE TWO p q taumax
  T=ONE[0]
  tau=[float(t)*1000/Fs for t in T] #convert to ms
  TMAX=taumax*1000/Fs #convert to ms
  pb.figure(1)
  pb.clf()
  pONE1,=pb.plot(tau,ONE[1],col1,linestyle='None',marker='.')
  pONE2,=pb.plot(tau,ONE[2],col2,linestyle='None',marker='.')
  pb.plot(tau,singlex(p,T,q,PoALL)[0][1],col1,tau,singlex(p,T,q,PoALL)[0][2],col2)
  pb.legend([pONE1, pONE2], ["$P_{one}$", "$P_{two}$"],fontsize=18)
  pb.grid(True)
  pb.xlim([0,TMAX])
  pb.ylim([0,1])
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlabel('t (ms)',fontsize=18)
  pb.ylabel('P(t)',fontsize=18)
  if directory !=None: pb.savefig(directory+'ONE.png')
  
  pb.figure(2)
  pb.clf()
  pb.grid(True)
  pTWO1,=pb.plot(tau,TWO[1],col1,linestyle='None',marker='.')
  pTWO2,=pb.plot(tau,TWO[2],col2,linestyle='None',marker='.')
  pb.plot(tau,singlex(p,T,q,PoALL)[1][1],col1,tau,singlex(p,T,q,PoALL)[1][2],col2)
  pb.legend([pTWO1, pTWO2], ["$P_{-}$", "$P_{o}$"],fontsize=18)
  pb.xlim([0,TMAX])
  pb.ylim([0,1])
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlabel('t (ms)',fontsize=18)
  pb.ylabel('P(t)',fontsize=18)
  if directory !=None: pb.savefig(directory+'TWO.png')
  
  if directory!=None: return
  else:
    pb.show()
  return

def get_scanparams(filename):
  x=np.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
  x=[i*10**9 for i in x]
  BIASstring=filename[-13]+'.'+filename[-10:-8]
  BIAS=float(BIASstring)
  NUMstring=filename[-7:-4]
  NUM=int(NUMstring)
  DIST=x[NUM]
  return (BIAS,DIST)

#-#-# double gaussian fit:
#-#-#

def twoGauss(p,x):
  A1, mu1, A2, mu2 = p
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  return G1+G2

def twoGaussErr(p, x, n): 
  error=twoGauss(p, x) - n
  if min(p[0],p[2])<0:
    error=1000*error
  #minsep=1.5
  #if np.abs(p[2]-p[5])<minsep or np.abs(p[8]-p[5])<minsep or np.abs(p[8]-p[3])<minsep:
  #  error=1000*error
  return error

def twoGaussFit(x,n,guess=(3000,20, 3000,25),mode=None):
  params,success=leastsq(twoGaussErr,guess,args=(x,n))
  A1,mu1,A2,mu2=params
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  if mode=='vocal':
    print 'Parameters for triple Gaussian fit to whole histogram:'
    print 'A1: '+str('%d'%A1)
    print 'sig1: '+str('%.2f'%sig1)
    print 'mu1: '+str('%.2f'%mu1)+'\n'
    print 'A2: '+str('%d'%A2)
    print 'sig2: '+str('%.2f'%sig2)
    print 'mu2: '+str('%.2f'%mu2)+'\n'
  
  tot=(params[0]+params[2])
  #make sure means are in the right order:
  if params[1]>params[3]:
    A1,sig1,mu1,A2,sig2,mu2=(A2,sig2,mu2,A1,sig1,mu1)
    params=(A1,mu1,A2,mu2)
  return params

#-#-# constrained 2-gauss fit:
#-#-#
def twoGauss_fixed(p, x, q):
  A1,A2=p
  mu1,mu2=q
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  return G1+G2

def twoGaussErr_fixed(p, x, n, q): 
  error=twoGauss_fixed(p, x, q) - n
  if p[0]<0 or p[1]<0: #if any amplitudes get negative
    error=1000*error
  #minsep=3
  #if np.abs(p[2]-p[5])<minsep or np.abs(p[8]-p[5])<minsep or np.abs(p[8]-p[3])<minsep: #if peaks get too close together
  #  error=1000*error
  return error

def twoGaussFit_fixed(x,n,q):
  p0=(3000, 3000)
  params,success=leastsq(twoGaussErr_fixed,p0,args=(x,n,q))
  return params
  
#-#-# single exponential fit: (used when there are two states)
#-#-#
def singlex(p,T,q,PoALL):#p=(e,) ; T is the array of tau values ; q is the additional parameter as a tuple, q=(gamma,)
  
  e,=p
  gamma,=q
  f=gamma*e
  
  M=np.array([[-e,f],[e,-f]])
  J,S=np.linalg.eig(M)
  
  ### one
  Q0=np.linalg.solve(S,PoALL[0])
  ### --- project negative state forward in time
  Pone=[];Ptwo=[]
  for t in T:
    Pt=np.dot(S,np.exp(t*J)*Q0)
    Pone.append(Pt[0]) #one
    Ptwo.append(Pt[1]) #two
  ONE=(T,Pone,Ptwo)

  ### two
  Q0=np.linalg.solve(S,PoALL[1])
  ### --- project negative state forward in time
  Pone=[];Ptwo=[]
  for t in T:
    Pt=np.dot(S,np.exp(t*J)*Q0)
    Pone.append(Pt[0]) #one
    Ptwo.append(Pt[1]) #two
  TWO=(T,Pone,Ptwo)
  return (ONE,TWO)

def singlexErr(p,T,dat,q,PoALL):
  ONE,TWO=singlex(p,T,q,PoALL)
  y=ONE[1]+ONE[2]+TWO[1]+TWO[2]
  return np.array(y)-dat

def singlexFit(T,dat,q,PoALL):
  p0=(0.376,)#(e0,)
  params,success=leastsq(singlexErr,p0,args=(T,dat,q,PoALL))
  return params

