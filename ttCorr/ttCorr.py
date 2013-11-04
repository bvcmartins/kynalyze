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

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

numbin=300
Imax=100
Fs=1e4

class TraceSet:
  def __init__(self):
        self.Pneg = np.nan
        self.Pneu = np.nan
        self.Ppos = np.nan
        self.a = np.nan
        self.b = np.nan
        self.c = np.nan
        self.d = np.nan
        self.Aneg = np.nan
        self.SIGneg = np.nan
        self.MUneg = np.nan
        self.Aneu = np.nan
        self.SIGneu = np.nan
        self.MUneu = np.nan
        self.Apos = np.nan
        self.SIGpos = np.nan
        self.MUpos = np.nan
        self.V = np.nan
        self.D = np.nan
        self.filname=''

def sigofmu(mu):
  sig=0.5+(1./20)*mu
  return sig

# get data
def getItrace(filename):
  f=open(filename,'rU')
  cond=True
  for x in range(35):f.readline()   #skip header
  line='asdfasdf'
  data=f.readlines() 
  Ifwd=[dt.split()[1] for dt in data]
  Ibwd=[dt.split()[2] for dt in data]
  Ibwd.reverse()
  I=[]
  I.extend(Ifwd)
  I.extend(Ibwd)
  I=[float(i)*1E12 for i in I]    #change units to pA
  return I

# clean data
def cleanFT(I,highFreq=1.5):
  # Current Trace:
  Fs=1e4
  T=float(len(I))/Fs
  t=np.linspace(0,T,len(I))
  # Fourier Transform:
  FT=scipy.fftpack.fft(np.array(I))#,100)
  freq=np.arange(float(len(FT)/2+1))/T
  #highFreq=1.5
  highInd=min(range(len(freq)), key=lambda i: abs(freq[i]-highFreq))
  FTkeep=[]
  FTkeep.extend(FT)
  for x in range(highInd)[1:]:
    FTkeep[x]=0
    FTkeep[-x]=0
  FTcut=[]
  FTcut.extend(FT)
  FTcut[0]=0
  for x in range(len(FT))[highInd:-highInd]:
    FTcut[x]=0
  
  Ikeep=scipy.fftpack.ifft(FTkeep)
  Icut=scipy.fftpack.ifft(FTcut)
  return Ikeep

#make histogram data
def getHist(I):
  n,bins=np.histogram(I,bins=np.linspace(0,Imax,numbin+1))
  ret=(n,bins)
  return ret

#show trace
def showTrace(I):
  t=[float(i)*1000/Fs for i in range(len(I))]
  pb.clf()
  pb.plot(t,I,'k',marker='.',linewidth=0.4)
  pb.grid(True)
  pb.xlim([200,400])
  pb.ylim([0,40])
  pb.xlabel('t (ms)')
  pb.ylabel('I (pA)')
  pb.show()
  return

#show histogram
def showHist(x,n,params):
  pb.clf()
  pb.bar(x,n,float(Imax)/numbin,color='0.8',linewidth=0.4,align='center')
  pb.plot(np.linspace(0,100,500),oneGauss(params[0:2],np.linspace(0,100,500)),color=col1,linewidth=1)
  pb.plot(np.linspace(0,100,500),oneGauss(params[2:4],np.linspace(0,100,500)),color=col2,linewidth=1)
  pb.plot(np.linspace(0,100,500),oneGauss(params[4:6],np.linspace(0,100,500)),color=col3,linewidth=1)
  pb.plot(np.linspace(0,100,500),threeGauss(params,np.linspace(0,100,500)),'--k',linewidth=1)
  pb.xlabel('I (pA)',fontsize=18)
  pb.ylabel('Counts',fontsize=18)
  pb.xlim([0,60])
  pb.grid(True)
  pb.show()
  return

# get time-correlation hist
def corrHist(I,Irange,tau):
  lower=Irange[0]
  upper=Irange[1]
  Icorr=[]
  k=tau
  for i in I[tau:]:
    if ((I[k-tau]>lower)and(I[k-tau]<upper)):
      Icorr.append(i)
    k+=1
  if len(Icorr)==0:
    return # can't do nothing with nothing
  n,bins=np.histogram(Icorr,bins=np.linspace(0,Imax,numbin+1))
  ret=(n,bins)
  return ret

def saveCorrHist(I,Irange,tau):
  n,bins=getHist(I)
  x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
  pb.clf()
  pb.bar(x,n,float(Imax)/numbin,color='0.8',linewidth=0.4,align='center')
  pb.plot(np.linspace(0,40,500),oneGauss(params[0:2],np.linspace(0,40,500)),color=col1,linewidth=1)
  pb.plot(np.linspace(0,40,500),oneGauss(params[2:4],np.linspace(0,40,500)),color=col2,linewidth=1)
  pb.plot(np.linspace(0,40,500),oneGauss(params[4:6],np.linspace(0,40,500)),color=col3,linewidth=1)
  pb.plot(np.linspace(0,40,500),threeGauss(params,np.linspace(0,40,500)),'--k',linewidth=1)
  pb.xlim([0,40])
  pb.grid(True)
  pb.show()
  return
  
#show correlation trace
def showTraceCorr(I,Irange,tau):
  t=[float(i)*1000/Fs for i in range(len(I))]
  lower=Irange[0]
  upper=Irange[1]
  Icorr=[]
  tcorr=[]
  k=tau
  for i in I[tau:]:
    if ((I[k-tau]>lower)and(I[k-tau]<upper)):
      Icorr.append(i)
      tcorr.append(t[k])
    k+=1
  if len(Icorr)==0:
    return # can't do nothing with nothing
  pb.clf()
  pb.plot(t,I,'k',marker='.',linewidth=0.4)
  pb.plot(tcorr,Icorr,col1,marker='.',linewidth=0)
  pb.grid(True)
  pb.xlim([200,400])
  pb.ylim([0,40])
  pb.xlabel('t (ms)')
  pb.ylabel('I (pA)')
  pb.text(360,37,'tau=%.1f'%(tau*1000/Fs))
  pb.show()
  return

#show correlation trace
def showTraceCorrs(I,Irange,taumax):
  t=[float(i)*1000/Fs for i in range(len(I))]
  lower=Irange[0]
  upper=Irange[1]
  for tau in range(taumax):    
    Icorr=[]
    tcorr=[]
    k=tau
    for i in I[tau:]:
      if ((I[k-tau]>lower)and(I[k-tau]<upper)):
        Icorr.append(i)
        tcorr.append(t[k])
      k+=1
    if len(Icorr)==0:
      return # can't do nothing with nothing
    pb.clf()
    pb.plot(t,I,'k',marker='.',linewidth=0.4)
    pb.plot(tcorr,Icorr,col1,marker='.',linewidth=0)
    pb.grid(True)
    pb.xlim([200,400])
    pb.ylim([0,40])
    pb.xlabel('t (ms)')
    pb.ylabel('I (pA)')
    pb.text(360,37,'tau=%.1f'%(tau*1000/Fs))
    pb.savefig('./traceCorr/'+str(tau)+'.png')

def showCorrHists(I,Irange,taumax,params):
  q=(params[1],params[3],params[5])
  n,bins=getHist(I)
  x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
  params=threeGaussFit(x,n)
  T,P1,P2,P3=([],[],[],[])
  for tau in range(taumax+1)[1:]:
    print '\r'+str('%d%%'%(100.*tau/taumax)),
    sys.stdout.flush()
    nCorr,binsCorr=corrHist(I,Irange,tau)
    params2=threeGaussFit_fixed(x,nCorr,q)
    pb.clf()
    pb.bar(x,n,float(Imax)/numbin,color='0.8',linewidth=0.4,align='center')
    pb.bar(x,nCorr,float(Imax)/numbin,color=col1,linewidth=0.4,align='center')
    pb.plot(np.linspace(0,40,500),oneGauss((params2[0]+params[1]),np.linspace(0,40,500)),color=col1,linewidth=1)
    pb.plot(np.linspace(0,40,500),oneGauss((params2[1]+params[3]),np.linspace(0,40,500)),color=col2,linewidth=1)
    pb.plot(np.linspace(0,40,500),oneGauss((params2[2]+params[5]),np.linspace(0,40,500)),color=col3,linewidth=1)
    #pb.plot(np.linspace(0,40,500),threeGauss(params,np.linspace(0,40,500)),'--k',linewidth=1)  
    pb.savefig('./corrHists/'+str('%3d'%tau)+'.png')
  return

# NOT USED FOR NOW (make a mode for this)
def PLOTS_corrHist(dirname,I,Irange,taumax):
  #print Irange
  if dirname[-1]!='/':
    dirname=dirname+'/'
  for tau in range(taumax+1):
    print 'working on tau='+str('%d'%tau)
    fname=dirname+str('%d'%Irange[0])+'pA'+'to'+str('%d'%Irange[1])+'pA'+str('-%03d'%tau)+'.png'
    n,bins=getHist(I)
    pb.clf()
    pb.bar(bins[1:],n,float(Imax)/numbin,color='0.5',alpha=0.5,linewidth=0)
    nCorr,binsCorr=corrHist(I,Irange,tau)
    pb.bar(binsCorr[1:],nCorr,float(Imax)/numbin,color=col3,linewidth=0.4)
    pb.grid(True)
    pb.xlabel('I (pA)',fontsize=16)
    pb.xlabel('counts',fontsize=16)
    plt.text(33,0.9*max(n),'tau='+str('%d'%tau))
    pb.savefig(fname)
  return

# NOT USED FOR NOW (make a mode for this) 
def threePLOTS_corrHist(dirname,I,Irange1,Irange2,Irange3,taumax):
  if dirname[-1]!='/':
    dirname=dirname+'/'
  for tau in range(taumax+1):
    print 'working on tau='+str('%d'%tau)
    fname=dirname+str('%03d'%tau)+'.png'
    n,bins=getHist(I)
    pb.clf()
    pb.bar(bins[1:],n,float(Imax)/numbin,color='0.5',alpha=0.5,linewidth=0)
    nCorr1,binsCorr1=corrHist(I,Irange1,tau)
    nCorr2,binsCorr2=corrHist(I,Irange2,tau)
    nCorr3,binsCorr3=corrHist(I,Irange3,tau)
    pb.bar(binsCorr1[1:],nCorr1,float(Imax)/numbin,color=col1,linewidth=0.4,alpha=0.5)
    pb.bar(binsCorr2[1:],nCorr2,float(Imax)/numbin,color=col2,linewidth=0.4,alpha=0.5)
    pb.bar(binsCorr3[1:],nCorr3,float(Imax)/numbin,color=col3,linewidth=0.4,alpha=0.5)
    pb.grid(True)
    pb.xlabel('I (pA)',fontsize=16)
    pb.xlabel('counts',fontsize=16)
    plt.text(33,0.9*max(n),'tau='+str('%d'%tau))
    pb.savefig(fname)
  return

def get_scanparams(filename):
  x=np.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
  x=[i*10**9 for i in x]
  BIASstring=filename[-13]+'.'+filename[-10:-8]
  BIAS=float(BIASstring)
  NUMstring=filename[-7:-4]
  NUM=int(NUMstring)
  DIST=x[NUM-1]
  return (BIAS,DIST)
  
def get_guessparams(BIAS,DIST):
  Q=0.007
  gamma=1.7
  PoNEG=np.exp(-2*DIST*gamma)/(Q**2+Q*np.exp(-DIST*gamma)+np.exp(-2*DIST*gamma))
  PoNEU=np.exp(-DIST*gamma)/(Q+np.exp(-DIST*gamma)+np.exp(-2*DIST*gamma)/Q)
  PoPOS=1/(1+np.exp(-DIST*gamma)/Q+np.exp(-2*DIST*gamma)/Q**2)
  mu1=6+(15./4.464105)*DIST+3.*np.sin(2*np.pi*DIST/0.72-1.2*np.pi)
  mu2=21+(5./4.464105)*DIST+3.*np.sin(2*np.pi*DIST/0.72-1.2*np.pi)
  mu3=40-(10./4.464105)*DIST+3.*np.sin(2*np.pi*DIST/0.72-1.2*np.pi)
  print PoNEG
  print PoNEU
  print PoPOS
  guessparams=(PoNEG*2*Fs*Imax/numbin, mu1, PoNEU*2*Fs*Imax/numbin, mu2, PoPOS*2*Fs*Imax/numbin, mu3)
  return guessparams

#-#-# single gaussian fit:
#-#-#
def oneGauss(p,x):
  A, mu = p
  sig=sigofmu(mu)
  G=(A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig))
  return G

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

def twoGaussFit(x,n,guess=(3000, 20, 3000, 25),mode=None):
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

#-#-# triple gaussian fit:
#-#-#
def threeGauss(p,x):
  A1, mu1, A2, mu2, A3, mu3 = p
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  sig3=sigofmu(mu3)
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

def threeGaussFit(x,n,guess=(3000, 20, 3000, 25, 3000, 30),mode=None):
  # This handles plots 7 through 14 well, but not outside of that
  # --> need better initial guesses 
  params,success=leastsq(threeGaussErr,guess,args=(x,n))
  A1,mu1,A2,mu2,A3,mu3=params
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  sig3=sigofmu(mu3)
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
  
#-#-# constrained 3-gauss fit:
#-#-#
def threeGauss_fixed(p, x, q):
  A1,A2,A3=p
  mu1,mu2,mu3=q
  sig1=sigofmu(mu1)
  sig2=sigofmu(mu2)
  sig3=sigofmu(mu3)
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  G3=(A3/(sig3*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu3)*(x-mu3)/(sig3*sig3))
  return G1+G2+G3

def threeGaussErr_fixed(p, x, n, q): 
  error=threeGauss_fixed(p, x, q) - n
  if p[0]<0 or p[1]<0 or p[2]<0: #if any amplitudes get negative
    error=1000*error
  #minsep=3
  #if np.abs(p[2]-p[5])<minsep or np.abs(p[8]-p[5])<minsep or np.abs(p[8]-p[3])<minsep: #if peaks get too close together
  #  error=1000*error
  return error

def threeGaussFit_fixed(x,n,q):
  p0=(3000, 3000, 3000)
  params,success=leastsq(threeGaussErr_fixed,p0,args=(x,n,q))
  return params

#-#-#-#-#-#-# 
#-#-#-#-#-#-# 
def analyze(filename):
  # read file into a list, I
  I=getItrace(filename) # can add: I=cleanFT(I,highFreq=2)
    
  # get and fit histogram of data
  n,bins=getHist(I)
  x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
  params=threeGaussFit(x,n)
  tot=params[0]+params[3]+params[6]
  totalCounts=tot*numbin/Imax
  Pneg=params[0]/tot #negative steady state
  Pneu=params[3]/tot #neutral steady state
  Ppos=params[6]/tot #positive steady state
  
  print 'Steady state probabilities (from triple Gaussian fit):'
  print 'Pneg: '+str('%.4f'%Pneg)
  print 'Pneu: '+str('%.4f'%Pneu)
  print 'Ppos: '+str('%.4f'%Ppos)+'\n'
  
  alpha=Pneg/Pneu
  beta=Pneu/Ppos
  
  if np.abs(totalCounts-len(I))>(float(len(I))/100):
    MESSAGE='---> WARNING: totalCounts is ' +str(totalCounts)
    print '-'*len(MESSAGE)
    print MESSAGE
    print '-'*len(MESSAGE)+'\n'
  
  # get time correlations over three ranges
  taumax=200
  IrangeNEG,IrangeNEU,IrangePOS,PoALL=ranges(params)
  print 'Working on negative range:'
  NEG=getCorr(I,IrangeNEG,taumax,params) # NEG=(T,P1,P2,P3)
  print '\nWorking on neutral range:'
  NEU=getCorr(I,IrangeNEU,taumax,params) # NEU=(T,P1,P2,P3)
  print '\nWorking on positive range:'
  POS=getCorr(I,IrangePOS,taumax,params) # POS=(T,P1,P2,P3)
  print '\n'
  
  q=(alpha,beta)
  p=dublexFit(np.array(NEG[0]),np.array(NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]),q,PoALL)
  a,b=p
  a=a*Fs # converts to units of Hz
  b=b*Fs # converts to units of Hz
  c,d=(a*alpha,b*beta)
  return (a,b,c,d)

def main():
  MESSAGE="ttCorr.py is a library of functions that are to be called by Analyze.py and SetAnalyze.py"
  print '\n'+'-'*min(80,len(MESSAGE))
  print '*'*min(80,len(MESSAGE))
  print '-'*min(80,len(MESSAGE))+'\n'
  print MESSAGE
  print '\n'+'-'*min(80,len(MESSAGE))
  print '*'*min(80,len(MESSAGE))
  print '-'*min(80,len(MESSAGE))+'\n'

#-#-#-#-#-#-#

if __name__ == '__main__':
  main()

# TO DO:
# - return error if failed to fit well
# - provide good guesses
# - try double and single gauss fits when triple fails... maybe even start with these


