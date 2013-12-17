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
import threestate as three

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

class TRACE:
  '''
  This class is the basic structure that will hold the data for a current trace. We start
  by simply assigning it a filename. Through Analyze.py, we fill in the rest of the data.
  It is important that there not be any "big data" in here, so that we can pass the object
  around without difficulty.
  
  DO NOT REMOVE STRUCTURES FROM THIS CLASS, OR CHANGE THEM! ONLY ADD THEM 
  '''
  def __init__(self,fname,numbin=300,Imax=100,Fs=1e4):
    
    MESSAGE='Initiating trace object for: '+fname;
    print '-'*len(MESSAGE)
    print MESSAGE
    print '-'*len(MESSAGE)+'\n'
    
    #basic info
    self.filename=fname
    self.dataset=''
    self.numbin=numbin
    self.Imax=Imax
    self.I = np.nan
    self.Itraj = np.nan
    if Fs==None:
      succ=False
      while succ==False:
        IN=raw_input("What is the sample rate? (e.g. 1e4) : ")
        try:
          #float(IN)
          self.Fs=float(IN)
          succ=True
        except ValueError:
          if IN!='':print "Whakina turkey jive is this?"
    else:
      self.Fs=Fs
        
    #experimental parameters
    self.bias = np.nan
    self.POStip = np.nan
    self.POSdb = np.nan
    self.dist = np.nan
    #histogram
    self.bins = np.nan
    self.x = np.nan
    self.n = np.nan
    #Gaussian fit
    self.histtype = '3' # should only have certain allowed values: ['3','2neg','2pos','1neg','1pos']
    self.params = np.nan
    ###the following should be derived from the above, not independently assigned
    self.Pneg = np.nan
    self.Pneu = np.nan
    self.Ppos = np.nan
    self.Aneg = np.nan
    self.SIGneg = np.nan
    self.MUneg = np.nan
    self.Aneu = np.nan
    self.SIGneu = np.nan
    self.MUneu = np.nan
    self.Apos = np.nan
    self.SIGpos = np.nan
    self.MUpos = np.nan
    #rates
    self.a = np.nan
    self.b = np.nan
    self.c = np.nan
    self.d = np.nan

####

# THIS IS A KLUDGE!
def sigofmu(mu):
  sig=0.5+(1./20)*mu
  return sig
# THAT WAS A KLUDGE!

def get_scanparams(TR):
  if open(TR.filename,'rU').readline()[:3]=='Exp':
    x=np.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
    x=[i*10**9 for i in x]
    BIASstring=TR.filename[-13]+'.'+TR.filename[-10:-8]
    BIAS=float(BIASstring)
    NUMstring=TR.filename[-7:-4]
    NUM=int(NUMstring)
    DIST=x[NUM-1]
    TR.bias=BIAS
    TR.dist=DIST
  else:
    f=open(TR.filename[:-7]+'.dat','rU')
    NUMPT=float(TR.filename[-7:-4])
    START=f.readline().split('      ')[1:]
    START=[float(start) for start in START]
    END=f.readline().split('      ')[1:]
    END=[float(end) for end in END]
    NUMPTS=int(f.readline().split('  ')[1])
    BIAS=float(f.readline().split()[1])
    TR.POStip=[START[0]+(END[0]-START[0])*NUMPT/NUMPTS,START[1]+(END[1]-START[1])*NUMPT/NUMPTS]
    TR.POSdb=END
    TR.bias=BIAS
    TR.dist=np.sqrt(sum((np.array(TR.POStip)-np.array(TR.POSdb))**2))
  return TR

# get data
def getItrace(TR):
  f=open(TR.filename,'rU')
  cond=True
  if f.readline()[:3]=='Exp':
    for x in range(35):f.readline()   #skip header
    data=f.readlines() 
    Ifwd=[dt.split()[1] for dt in data]
    Ibwd=[dt.split()[2] for dt in data]
    Ibwd.reverse()
    I=[]
    I.extend(Ifwd)
    I.extend(Ibwd)
    I=[float(i)*1E12 for i in I]    #change units to pA
  else: #assume that it's one of my simple data files
    data=f.readlines()
    I=[float(dt.split()[1])*1E12 for dt in data]
  TR.I=I
  return TR

#make histogram data
def getHist(TR):
  n,bins=np.histogram(TR.I,bins=np.linspace(0,TR.Imax,TR.numbin+1))
  TR.n=n
  TR.bins=bins
  TR.x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
  return TR

def histFit(TR):
  guessparams=get_guessparams(TR)
  TR.params=three.threeGaussFit(TR,mode='vocal',guess=guessparams)
  showHist(TR)
  
  ### Fitting the histogram with Gaussians: 
  success=False
  while success==False:
    var=raw_input("OK or adjust guess ('y', or '3' to change initial guess) [add two state capability later] ?")
    if var=='y':
      success=True
      tot=TR.params[0]+TR.params[2]+TR.params[4]
      totalCounts=tot*TR.numbin/TR.Imax
      TR.Pneg=TR.params[0]/tot #negative steady state
      TR.Pneu=TR.params[2]/tot #neutral steady state
      TR.Ppos=TR.params[4]/tot #positive steady state
      TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=TR.params   
      TR.SIGneg,TR.SIGneu,TR.SIGpos=(sigofmu(TR.params[1]),sigofmu(TR.params[3]),sigofmu(TR.params[5]))
    elif var=='3':
      three.adjustGuess(TR)
      TR.params=three.threeGaussFit(TR,mode='vocal',guess=guessparams)
      showHist(TR)
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y':
          succ=True
          success=True
          TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=TR.params
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(sigofmu(TR.params[1]),sigofmu(TR.params[3]),sigofmu(TR.params[5]))
          tot=TR.params[0]+TR.params[2]+TR.params[4]
          totalCounts=tot*TR.numbin/TR.Imax
          TR.Pneg=TR.params[0]/tot #negative steady state
          TR.Pneu=TR.params[2]/tot #neutral steady state
          TR.Ppos=TR.params[4]/tot #positive steady state
        elif var=='n':
          succ=True
          success=False
        elif var!='':
          print "Ain't nobody got time for that!"
    elif var!='':
      print "give me something I can work with!"
  print 'Steady state probabilities (from multiple Gaussian fit):'
  print 'Pneg: '+str('%.4f'%TR.Pneg)
  print 'Pneu: '+str('%.4f'%TR.Pneu)
  print 'Ppos: '+str('%.4f'%TR.Ppos)+'\n'
  if np.abs(totalCounts-sum(TR.n))>(float(sum(TR.n))/100):
    MESSAGE='---> WARNING: totalCounts is ' +str(totalCounts)
    print '-'*len(MESSAGE)
    print MESSAGE
    print '-'*len(MESSAGE)+'\n'
  return TR

def get_guessparams(TR):
  #Q=0.007
  #gamma=1.7
  PoNEG=0.33
  PoNEU=0.33
  PoPOS=0.33
  mu2=TR.x[list(TR.n).index(max(TR.n))]
  mu1=0.9*mu2
  mu3=1.1*mu2  
  guessparams=(PoNEG*2*TR.Fs*TR.Imax/TR.numbin, mu1, PoNEU*2*TR.Fs*TR.Imax/TR.numbin, mu2, PoPOS*2*TR.Fs*TR.Imax/TR.numbin, mu3)
  return guessparams

#show histogram
def showHist(TR,pars=None): # leave these inputs (don't just put TR, because we also use it for guesses)
  if pars!=None: params=pars
  else: params=TR.params
  xs=[]
  for i in range(len(TR.x)):
    if float(TR.n[i])/max(TR.n) > 0.01:
      xs.append(TR.x[i])
  Xdiff=max(xs)-min(xs)
  Xav=0.5*(max(xs)+min(xs))
  Xmin=Xav-0.75*Xdiff
  Xmax=Xav+0.75*Xdiff
  pb.clf()
  pb.bar(TR.x,TR.n,float(TR.Imax)/TR.numbin,color='0.8',linewidth=0.4,align='center')
  pb.plot(np.linspace(0,100,500),oneGauss(params[0:2],np.linspace(0,100,500)),color=col1,linewidth=1)
  pb.plot(np.linspace(0,100,500),oneGauss(params[2:4],np.linspace(0,100,500)),color=col2,linewidth=1)
  pb.plot(np.linspace(0,100,500),oneGauss(params[4:6],np.linspace(0,100,500)),color=col3,linewidth=1)
  pb.plot(np.linspace(0,100,500),three.threeGauss(params,np.linspace(0,100,500)),'--k',linewidth=1)
  pb.xlabel('I (pA)',fontsize=20)
  pb.ylabel('Counts',fontsize=20)
  pb.xticks(fontsize=16)
  pb.yticks(fontsize=16)
  pb.xlim([Xmin,Xmax])
  pb.grid(True)
  pb.show()
  return

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

def getItraj(TR,taumax):
  taumax=int(taumax*TR.Fs+1)
  print taumax
  Itraj=[]
  for k in range(len(TR.I)-taumax):
    Itraj.append(TR.I[k:k+taumax])
  Itraj.sort(key=lambda traj: traj[0])
  TR.Itraj=Itraj
  return

#-#-# single gaussian fit:
#-#-#
def oneGauss(p,x):
  A, mu = p
  sig=sigofmu(mu)
  G=(A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig))
  return G

#-#-#-#-#-#-# 
#-#-#-#-#-#-# 

def Analyze(filename):
  TR=TRACE(filename) # initialize the trace object
  TR=get_scanparams(TR) # get bias and pos and dist
  TR=getItrace(TR) # add the data to TR.I
  TR=getHist(TR) # get histogram
  TR=histFit(TR) # fit the histogram
  
  taumax=10e-3 # 10 ms
  # get time correlations over three ranges
  getItraj(TR,taumax)
  IrangeNEG,IrangeNEU,IrangePOS,PoALL=three.ranges(TR.params)
  print 'Negative range:'
  print IrangeNEG
  print 'Neutral range:'
  print IrangeNEU
  print 'Positive range:'
  print IrangePOS
  alpha=TR.Pneg/TR.Pneu
  beta=TR.Pneu/TR.Ppos
  q=(alpha,beta)
  print '\nWorking on negative range:'
  NEG=three.getCorr(TR,IrangeNEG,taumax) # NEG=(T,P1,P2,P3)
  print '\nWorking on neutral range:'
  NEU=three.getCorr(TR,IrangeNEU,taumax) # NEU=(T,P1,P2,P3)
  print '\nWorking on positive range:'
  POS=three.getCorr(TR,IrangePOS,taumax) # POS=(T,P1,P2,P3)
  print '\n'
  
  p=three.dublexFit(np.array(NEG[0]),np.array(NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]),q,PoALL)
  a,b=p
  TR.a=a#*TR.Fs # converts to units of Hz --- no longer necessary
  TR.b=b#*TR.Fs # converts to units of Hz --- no longer necessary
  TR.c,TR.d=(a*alpha,b*beta)
  print 'a,b,c,d are:'
  print str('%.4f, '%TR.a)+str('%.4f, '%TR.b)+str('%.4f, '%TR.c)+str('%.4f, '%TR.d)
  #three.dublexPLOT(NEG,NEU,POS,p,q,PoALL,taumax)
  three.dublexPLOT(NEG,NEU,POS,(500,500),q,PoALL,taumax)
  return TR

def main():
  MESSAGE="rr.py is a library of functions that are to be called by Analyze.py and SetAnalyze.py"
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
