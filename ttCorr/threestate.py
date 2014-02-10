import sys
import os
import numpy as np
import pylab as pb
import scipy.fftpack
import scipy.special
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import ttCorr as ttC

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

def adjustGuess(TR):
  x=TR.x
  n=TR.n
  guessparams=TR.params
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
    elif 'A3=' in var: 
      guessparams=list(guessparams)
      guessparams[4]=float(var[3:])
      guessparams=tuple(guessparams)
    elif 'mu1=' in var: 
      guessparams=list(guessparams)
      guessparams[1]=float(var[4:])
      guessparams=tuple(guessparams)
    elif 'mu2=' in var: 
      guessparams=list(guessparams)
      guessparams[3]=float(var[4:])
      guessparams=tuple(guessparams)
    elif 'mu3=' in var: 
      guessparams=list(guessparams)
      guessparams[5]=float(var[4:])
      guessparams=tuple(guessparams)
    elif var=='show':
      ttC.showHist(TR,guessparams)
    elif var=='y':
      good=True
    elif var!='':
      print "Use one of the accepted parameters (A1,mu1,A2,mu2,A3,mu3)"
  TR.params=guessparams
  return None

def getCorr(TR,Irange,taumax):
  P=[]
  T,P1,P2,P3=([],[],[],[])
  for I in TR.x:
    p1=ttC.oneGauss((TR.params[0],TR.params[1]),I)
    p2=ttC.oneGauss((TR.params[2],TR.params[3]),I)
    p3=ttC.oneGauss((TR.params[4],TR.params[5]),I)
    p=np.array([p1,p2,p3])
    p=p/p.sum()
    P.append(p)
  if TR.Fs>1e4:
    TAU=range(0,int(taumax*TR.Fs)+1,int(TR.Fs/1e4))
  else:
    TAU=range(int(taumax*TR.Fs)+1)#[3:] # change 3 to 0 later
  lower=Irange[0]
  upper=Irange[1]
  INDlower=min(range(len(TR.I))[:-int(taumax*TR.Fs)-1], key=lambda k: np.abs(TR.Itraj[k][0]-lower))
  INDupper=min(range(len(TR.I))[:-int(taumax*TR.Fs)-1], key=lambda k: np.abs(TR.Itraj[k][0]-upper))
  Itrbl=np.array(TR.Itraj[INDlower:INDupper]) # this is a "trajectory block"
  for tau in TAU:
    print '\r'+str('%d%%'%(100.*tau/taumax/TR.Fs)),
    sys.stdout.flush()
    nCorr,binsCorr=np.histogram(Itrbl[:,tau],bins=np.linspace(0,TR.Imax,TR.numbin+1))
    ptot=np.zeros(3)
    ktot=0
    for k in range(len(TR.x)):
      ptot+=nCorr[k]*P[k]
      ktot+=nCorr[k]
    ptot=ptot/ktot
    T.append(float(tau)/TR.Fs)
    P1.append(ptot[0])
    P2.append(ptot[1])
    P3.append(ptot[2])
  return (T,P1,P2,P3)


def ranges(params):
  A1,mu1,A2,mu2,A3,mu3=params
  sig1=ttC.sigofmu(mu1)
  sig2=ttC.sigofmu(mu2)
  sig3=ttC.sigofmu(mu3)
  IrangeNEG=(0,mu2-2*sig2)
  IrangeNEU=(mu2-sig2,mu2+sig2)
  IrangePOS=(mu2+2*sig2,mu3+8*sig3)
  #NEG integrals
  NEGneg=A1*(scipy.special.erf((IrangeNEG[1]-mu1)/(np.sqrt(2)*sig1))-scipy.special.erf((IrangeNEG[0]-mu1)/(np.sqrt(2)*sig1)))
  NEGneu=A2*(scipy.special.erf((IrangeNEG[1]-mu2)/(np.sqrt(2)*sig2))-scipy.special.erf((IrangeNEG[0]-mu2)/(np.sqrt(2)*sig2)))
  NEGpos=A3*(scipy.special.erf((IrangeNEG[1]-mu3)/(np.sqrt(2)*sig3))-scipy.special.erf((IrangeNEG[0]-mu3)/(np.sqrt(2)*sig3)))
  tot=NEGneg+NEGneu+NEGpos
  PoNEG=[NEGneg/tot,NEGneu/tot,NEGpos/tot]
  #NEU integrals
  NEUneg=A1*(scipy.special.erf((IrangeNEU[1]-mu1)/(np.sqrt(2)*sig1))-scipy.special.erf((IrangeNEU[0]-mu1)/(np.sqrt(2)*sig1)))
  NEUneu=A2*(scipy.special.erf((IrangeNEU[1]-mu2)/(np.sqrt(2)*sig2))-scipy.special.erf((IrangeNEU[0]-mu2)/(np.sqrt(2)*sig2)))
  NEUpos=A3*(scipy.special.erf((IrangeNEU[1]-mu3)/(np.sqrt(2)*sig3))-scipy.special.erf((IrangeNEU[0]-mu3)/(np.sqrt(2)*sig3)))
  tot=NEUneg+NEUneu+NEUpos
  PoNEU=[NEUneg/tot,NEUneu/tot,NEUpos/tot]
  #POS integrals
  POSneg=A1*(scipy.special.erf((IrangePOS[1]-mu1)/(np.sqrt(2)*sig1))-scipy.special.erf((IrangePOS[0]-mu1)/(np.sqrt(2)*sig1)))
  POSneu=A2*(scipy.special.erf((IrangePOS[1]-mu2)/(np.sqrt(2)*sig2))-scipy.special.erf((IrangePOS[0]-mu2)/(np.sqrt(2)*sig2)))
  POSpos=A3*(scipy.special.erf((IrangePOS[1]-mu3)/(np.sqrt(2)*sig3))-scipy.special.erf((IrangePOS[0]-mu3)/(np.sqrt(2)*sig3)))
  tot=POSneg+POSneu+POSpos
  PoPOS=[POSneg/tot,POSneu/tot,POSpos/tot]
  
  PoALL=(PoNEG,PoNEU,PoPOS)
  return (IrangeNEG,IrangeNEU,IrangePOS, PoALL)

def dublexPLOT(NEG,NEU,POS,p,q,PoALL,taumax,directory=None):
#NEG NEU POS p q taumax
  tau=NEG[0]
  
  pb.figure(1)
  pb.clf()
  pNEG1,=pb.plot(np.array(tau)*1000,NEG[1],col1,linestyle='None',marker='.')
  pNEG2,=pb.plot(np.array(tau)*1000,NEG[2],col2,linestyle='None',marker='.')
  pNEG3,=pb.plot(np.array(tau)*1000,NEG[3],col3,linestyle='None',marker='.')
  pb.plot(np.array(tau)*1000,dublex(p,tau,q,PoALL)[0][1],col1,np.array(tau)*1000,dublex(p,tau,q,PoALL)[0][2],col2,np.array(tau)*1000,dublex(p,tau,q,PoALL)[0][3],col3)
  pb.legend([pNEG1, pNEG2, pNEG3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=18)
  pb.grid(True)
  pb.xlim([0,taumax*1000])
  pb.ylim([0,1])
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlabel('t (ms)',fontsize=18)
  pb.ylabel('P(t)',fontsize=18)
  if directory !=None: pb.savefig(directory+'NEG.png')
  
  pb.figure(2)
  pb.clf()
  pb.grid(True)
  pNEU1,=pb.plot(np.array(tau)*1000,NEU[1],col1,linestyle='None',marker='.')
  pNEU2,=pb.plot(np.array(tau)*1000,NEU[2],col2,linestyle='None',marker='.')
  pNEU3,=pb.plot(np.array(tau)*1000,NEU[3],col3,linestyle='None',marker='.')
  pb.plot(np.array(tau)*1000,dublex(p,tau,q,PoALL)[1][1],col1,np.array(tau)*1000,dublex(p,tau,q,PoALL)[1][2],col2,np.array(tau)*1000,dublex(p,tau,q,PoALL)[1][3],col3)
  pb.legend([pNEU1, pNEU2, pNEU3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=18)
  pb.xlim([0,taumax*1000])
  pb.ylim([0,1])
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlabel('t (ms)',fontsize=18)
  pb.ylabel('P(t)',fontsize=18)
  if directory !=None: pb.savefig(directory+'NEU.png')
  
  pb.figure(3)
  pb.clf
  pPOS1,=pb.plot(np.array(tau)*1000,POS[1],col1,linestyle='None',marker='.')
  pPOS2,=pb.plot(np.array(tau)*1000,POS[2],col2,linestyle='None',marker='.')
  pPOS3,=pb.plot(np.array(tau)*1000,POS[3],col3,linestyle='None',marker='.')
  pb.plot(np.array(tau)*1000,dublex(p,tau,q,PoALL)[2][1],col1,np.array(tau)*1000,dublex(p,tau,q,PoALL)[2][2],col2,np.array(tau)*1000,dublex(p,tau,q,PoALL)[2][3],col3)
  pb.legend([pPOS1, pPOS2, pPOS3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=18)
  pb.grid(True)
  pb.xlim([0,taumax*1000])
  pb.ylim([0,1])
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlabel('t (ms)',fontsize=18)
  pb.ylabel('P(t)',fontsize=18)
  if directory !=None: pb.savefig(directory+'POS.png')
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

#-#-# triple gaussian fit:
#-#-#
def threeGauss(p,x):
  A1, mu1, A2, mu2, A3, mu3 = p
  sig1=ttC.sigofmu(mu1)
  sig2=ttC.sigofmu(mu2)
  sig3=ttC.sigofmu(mu3)
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  G3=(A3/(sig3*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu3)*(x-mu3)/(sig3*sig3))
  return G1+G2+G3

def threeGaussErr(p, x, n): 
  error=threeGauss(p, x) - n
  if min(p[0],p[2],p[4])<0:
    error=1000*error
  #minsep=1.5
  #if np.abs(p[2]-p[5])<minsep or np.abs(p[8]-p[5])<minsep or np.abs(p[8]-p[3])<minsep:
  #  error=1000*error
  return error

def threeGaussFit(TR,guess=(3000, 20, 3000, 25, 3000, 30),mode=None):
  x=TR.x
  n=TR.n
  params,success=leastsq(threeGaussErr,guess,args=(x,n))
  A1, mu1, A2, mu2, A3, mu3 = params
  sig1=ttC.sigofmu(mu1)
  sig2=ttC.sigofmu(mu2)
  sig3=ttC.sigofmu(mu3)
  if mode=='vocal':
    print 'Parameters for triple Gaussian fit to whole histogram:'
    print 'A1: '+str('%d'%A1)
    print 'sig1: '+str('%.2f'%sig1)
    print 'mu1: '+str('%.2f'%mu1)+'\n'
    print 'A2: '+str('%d'%A2)
    print 'sig2: '+str('%.2f'%sig2)
    print 'mu2: '+str('%.2f'%mu2)+'\n'
    print 'A3: '+str('%d'%A3)
    print 'sig3: '+str('%.2f'%sig3)
    print 'mu3: '+str('%.2f'%mu3)+'\n'
  
  tot=(params[0]+params[2]+params[4])
  #make sure means are in the right order:
  if params[1]>params[3]:
    A1,sig1,mu1,A2,sig2,mu2=(A2,sig2,mu2,A1,sig1,mu1)
    params=(A1,mu1,A2,mu2,A3,mu3)
  if params[1]>params[5]:
    A1,sig1,mu1,A2,sig2,mu2,A3,sig3,mu3=(A3,sig3,mu3,A1,sig1,mu1,A2,sig2,mu2)
    params=(A1,mu1,A2,mu2,A3,mu3)
  elif params[3]>params[5]:
    A2,sig2,mu2,A3,sig3,mu3=(A3,sig3,mu3,A2,sig2,mu2)
    params=(A1,mu1,A2,mu2,A3,mu3)
  if params[1]>params[3] or params[1]>params[5] or params[3]>params[5]:
    print 'the means are not in order: you need to alter the code to get this right'
  return params

#-#-# double exponential fit: (used when there are three states)
#-#-#
def dublex(p,T,q,PoALL):#p=(a,b) ; T is the array of tau values ; q is the additional parameters as a tuple, q=(alpha,beta,...) #(T,a,b):
  a,b=p
  alpha,beta=q
  c=a*alpha
  d=b*beta
  
  M=np.array([[-a,c,0],[a,-b-c,d],[0,b,-d]])
  J,S=np.linalg.eig(M)
  #J=J[np.newaxis]
  #J=np.dot(J.T,np.array([[1,1,1]]))*np.identity(3)
  ### neg
  Q0=np.linalg.solve(S,PoALL[0])
  ### --- project negative state forward in time
  Pneg=[];Pneu=[];Ppos=[]
  for t in T:
    Pt=np.dot(S,np.exp(t*J)*Q0)
    #Pt=np.array(PoALL[0]*S*J*S.I)
    #Pt=Pt.flatten()
    Pneg.append(Pt[0]) #neg
    Pneu.append(Pt[1]) #neu
    Ppos.append(Pt[2]) #pos
  NEG=(T,Pneg,Pneu,Ppos)

  ### neu
  Q0=np.linalg.solve(S,PoALL[1])
  ### --- project negative state forward in time
  Pneg=[];Pneu=[];Ppos=[]
  for t in T:
    Pt=np.dot(S,np.exp(t*J)*Q0)
    Pneg.append(Pt[0]) #neg
    Pneu.append(Pt[1]) #neu
    Ppos.append(Pt[2]) #pos
  NEU=(T,Pneg,Pneu,Ppos)
  
  ### pos
  Q0=np.linalg.solve(S,PoALL[2])
  ### --- project negative state forward in time
  Pneg=[];Pneu=[];Ppos=[]
  for t in T:
    Pt=np.dot(S,np.exp(t*J)*Q0)
    Pneg.append(Pt[0]) #neg
    Pneu.append(Pt[1]) #neu
    Ppos.append(Pt[2]) #pos
  POS=(T,Pneg,Pneu,Ppos)
  return (NEG,NEU,POS)

def dublexErr(p,T,dat,q,PoALL):
  NEG,NEU,POS=dublex(p,T,q,PoALL)
  y=NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]
  return np.array(y)-dat

def dublexFit(T,dat,q,PoALL):
  p0=(500.,500.)#(a0,b0)
  #params,success=leastsq(dublexErr,p0,args=(T,dat,q,PoALL))
  def CVdublex(T,a,b):
    NEG,NEU,POS=dublex((a,b),T,q,PoALL)
    return NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]
  params,pcov=curve_fit(CVdublex,T,dat,p0=p0)
  errors=tuple([np.sqrt(pcov[k][k]) for k in range(len(params))])
  print "params are:"
  print params
  print "errors are:"
  print errors
  print "pcov is:"
  print pcov
  return params
  
