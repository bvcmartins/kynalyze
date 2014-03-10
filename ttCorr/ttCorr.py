import sys
import os
import copy
import numpy as np
import pylab as pb
import scipy.fftpack
import scipy.special
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
#import pickle
import matplotlib.pyplot as plt
import threestate as three
import tlib1 as t1
from operator import itemgetter
import re
import scipy.special

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'
colours=[col1,col2,col3,col4,'r','g','b',]

class TRACE:
  '''
  This class is the basic structure that will hold the data for a current trace. We start
  by simply assigning it a filename. Through Analyze.py, we fill in the rest of the data.
  It is important that there not be any "big data" in here, so that we can pass the object
  around without difficulty.
  '''
  def __init__(self,fname,numbin=300,Imax=100,Fs=1e4):
    
    MESSAGE='Initiating trace object for: '+fname;
    print '-'*len(MESSAGE)
    print MESSAGE
    print '-'*len(MESSAGE)+'\n'
    
    self.filename=fname
    
    ## directories
    motherdir=os.path.dirname(os.path.realpath(fname))
    kyndir=motherdir+'/kynalysis/'
    datdir=kyndir+'/kyndat/'
    figdir=kyndir+'/kynfig/'
    if not os.path.exists(kyndir): os.makedirs(kyndir)
    if not os.path.exists(datdir): os.makedirs(datdir)
    if not os.path.exists(figdir): os.makedirs(figdir)    
    self.kyndir=kyndir
    self.datdir=datdir
    self.figdir=figdir
    
    #basic info
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
    self.Nstates = None # should only have certain allowed values: ['3','2neg','2pos','1neg','1pos']
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
###### inserted by Bruno on 2013-12-20
    #smooth
    self.y = np.nan
    # peak finder
    self.peak = []
    self.amp = []
    self.popt = []
    self.sens = 10  # arbitrary parameter
    self.radius = []
    self.center = []
######

  def get_scanparams(self):
    if open(self.filename,'rU').readline()[:3]=='Exp':
      x=np.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
      x=[i*10**9 for i in x]
      BIASstring=self.filename[-13]+'.'+self.filename[-10:-8]
      BIAS=float(BIASstring)
      NUMstring=self.filename[-7:-4]
      NUM=int(NUMstring)
      DIST=x[NUM-1]
      self.bias=BIAS
      self.dist=DIST
    else:
      f=open(self.filename[:-7]+'.dat','rU')
      NUMPT=float(self.filename[-7:-4])
      START=f.readline().split('      ')[1:]
      START=[float(start) for start in START]
      END=f.readline().split('      ')[1:]
      END=[float(end) for end in END]
      NUMPTS=int(f.readline().split('  ')[1])
      BIAS=float(f.readline().split()[1])
      self.POStip=[START[0]+(END[0]-START[0])*NUMPT/NUMPTS,START[1]+(END[1]-START[1])*NUMPT/NUMPTS]
      self.POSdb=END
      self.bias=BIAS
      self.dist=np.sqrt(sum((np.array(self.POStip)-np.array(self.POSdb))**2))
    return None

  # get data
  def getItrace(self):
    f=open(self.filename,'rU')
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
#      self.I=[float(i)*1E12 for i in I] 
      self.I=[float(i) for i in I]    #change units to pA  - Bruno 14-03-06: removed 1E12.
    else: #assume that it's one of my simple data files
      data=f.readlines()
#      self.I=[float(dt.split()[1])*1E12 for dt in data]
      self.I=[float(dt.split()[1]) for dt in data]  # Bruno 14-03-06: removed 1E12.

    return None

  #make histogram data
  def getHist(self):
    n,bins=np.histogram(self.I,bins=np.linspace(min(self.I),max(self.I),0.3*self.numbin+1)) # Bruno 14-03-05: modified
    self.n=n
    self._n=copy.deepcopy(n) # this is a backup so that we can reset n later. It means we can smooth or alter n as much as we want and still return to the raw histogram.
    self.bins=bins
    self.x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
    return None
  
  #show histogram
  def plotHist(self,pars=None,SHOW=True,SAVE=True): # leave these inputs (don't just put TR, because we also use it for guesses)
    skipfits=False
    if pars!=None: params=pars
    else:
      if not np.isnan(np.sum(self.params)): params=self.params
      else: skipfits=True
    xs=[]
    for i in range(len(self.x)):
      if float(self.n[i])/max(self.n) > 0.01:
        xs.append(self.x[i])
    Xdiff=max(xs)-min(xs)
    Xav=0.5*(max(xs)+min(xs))
    Xmin=Xav-0.75*Xdiff
    Xmax=Xav+0.75*Xdiff
    pb.clf()
    pb.bar(self.x,self.n,float(self.Imax)/self.numbin,color='0.8',linewidth=0.4,align='center')
    if not skipfits:
      for i in range(len(params)/3):
        pb.plot(np.linspace(0,100,500),nGauss(np.linspace(0,100,500),params[i*3:i*3+3]),color=colours[i],linewidth=1)
      pb.plot(np.linspace(0,100,500),nGauss(np.linspace(0,100,500),params),'--k',linewidth=1)
    pb.xlabel('I (pA)',fontsize=20)
    pb.ylabel('Counts',fontsize=20)
    pb.xticks(fontsize=16)
    pb.yticks(fontsize=16)
    pb.xlim([Xmin,Xmax])
    pb.grid(True)
    if SAVE:
      fn=os.path.basename(self.filename)
      pb.savefig(self.figdir+'hist_'+fn[:-4]+'.eps')
      f=open(self.datdir+fn[:-4]+'.dat','w')
      f.write('#hist\n')
      f.write('BinCentre(pA)\tBinMin(pA)\tBinMax(pA)\tCounts\n')
      for i in range(len(self.n)):
        f.write(str(self.x[i])+'\t'+str(self.bins[i])+'\t'+str(self.bins[i+1])+'\t'+str(self.n[i])+'\n')
      f.close
    else: self.plotHist(SHOW=SHOW)
    #if savename!=None:pb.savefig(savename)
    if SHOW:pb.show()
    return
    
  def resetHist(self):
    self.n=copy.deepcopy(self._n)
    return

####

def sigofmu(mu):
  sig=0.5+(1./20)*mu
#  sig=0.05+(1./20)*mu
  return sig

def histFit(TR):
  guessparams=get_guessparams(TR)
  TR.params=three.threeGaussFit(TR,mode='vocal',guess=guessparams)
  #TR.plotHist()
  return # temporary return - inserted by Bruno on 2013-12-20
 
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
      TR.plotHist()
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
  
def nGauss(x,params):
  G=[]
  for i in range(len(params)/3):
    mu,A,sig=params[3*i:3*i+3]
    G.append((A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig)))
  return sum(G)

#-#-#-#-#-#-# FUNCTIONS FROM tlib1.py:
#-#-#-#-#-#-# 

def smooth(TR):
  window_len=11
  window='hanning'
  s=np.r_[TR.n[window_len-1:0:-1],TR.n,TR.n[-1:-window_len:-1]]
  w=eval('np.'+window+'(window_len)')
  TR.y=np.convolve(w/w.sum(),s,mode='valid')
  TR.y=TR.y[5:-5]

#  pb.bar(TR.x,TR.n,float(TR.Imax)/TR.numbin,color='0.8',linewidth=0.4,align='center')
  pb.plot(TR.x,TR.y,color=col1,linewidth=1)
#  pb.show() 

def peak_finder(TR):
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]
  for i in np.where(pk==1)[0]:
    if TR.y[int(i)]>TR.sens:
      TR.peak.append(TR.x[int(i)])
      TR.amp.append(TR.y[int(i)])

  for i in TR.amp:
    k=int(np.where(TR.y==i)[0])
    p=TR.x[k]
    c=TR.y[k]/np.sqrt(2.)
    while True:
      if TR.y[k]<c:
        break
      k+=1
    TR.radius.append(abs(p-TR.x[k]))

  print TR.peak
  print TR.amp
  print TR.radius

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
  return error

def oneGauss2(x,mu1=0.0,A1=0.0,sig1=0.0):
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  return G1

def twoGauss2(x,mu1=0.0,A1=0.0,sig1=0.0,mu2=0.0,A2=0.0,sig2=0.0):
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  return G1+G2

def threeGauss2(x,mu1=0.0,A1=0.0,sig1=0.0,mu2=0.0,A2=0.0,sig2=0.0,mu3=0.0,A3=0.0,sig3=0.0):
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  G3=(A3/(sig3*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu3)*(x-mu3)/(sig3*sig3))
  return G1+G2+G3

def nGauss2(x,params):
  G=[]
  for i in range(len(params)/3):
    mu,A,sig=params[3*i:3*i+3]
    G.append((A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig)))
  return sum(G)

def threeGaussFit(TR):
  print len(TR.peak)
 
  if len(TR.peak)>0:
    if len(TR.peak)==1:
      sig=[TR.radius[0]]
      popt,pcov=curve_fit(oneGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0]])
      TR.popt.append(popt[0])
      TR.popt.append(popt[1])
      TR.popt.append(popt[2])
#      pb.plot(TR.x,oneGauss2(TR.x,popt[0],popt[1],popt[2]),color=col4,linewidth=1)
#      pb.show()
      TR.n=TR.n-oneGauss2(TR.x,popt[0],popt[1],popt[2])

    elif len(TR.peak)==2:
      sig=[TR.radius[0],TR.radius[1]]
      popt,pcov=curve_fit(twoGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1]])
      for i in range(6):
        TR.popt.append(popt[i])
  
#      pb.plot(TR.x,twoGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),color=col4,linewidth=1)
#      pb.show()
      TR.n=TR.n-twoGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

    else:
      sig=[TR.radius[0],TR.radius[1],TR.radius[2]]
      popt,pcov=curve_fit(threeGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1],TR.peak[2],TR.amp[2],sig[2]])
      for i in range(9):
        TR.popt.append(popt[i])
  
      pb.plot(TR.x,threeGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8]),color=col4,linewidth=1)
      pb.show()
      TR.n=TR.n-threeGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8])

  TR.sens=0.02*max(TR.y)

def hidden_peak(TR,filename):
  window_len=11
  window='hanning'
  s=np.r_[TR.n[window_len-2:0:-1],TR.n,TR.n[-1:-window_len:-1]]
  w=eval('np.'+window+'(window_len)')
  TR.y=np.convolve(w/w.sum(),s,mode='valid')
  TR.y=TR.y[4:-5]
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]
 
  x=TR.x
  n=TR.n
 
  peak=[]
  amp=[]
  radius=[]
  radius2=[]

#  pb.plot(TR.x,TR.y,color=col1,linewidth=1)
#  pb.show() 


  for i in np.where(pk==1)[0]:
    min=TR.x[len(TR.x)-1]
    if TR.x[int(i)] != TR.x[0] and TR.x[int(i)] != TR.x[len(TR.x)-1] and TR.y[int(i)]>TR.sens:
      for j in TR.peak:
        if abs(TR.x[int(i)]-j) < min :
          min = abs(TR.x[int(i)]-j) 
          min2 = 2*TR.radius[TR.peak.index(j)]
      if min > min2:
        peak.append(TR.x[int(i)])
        amp.append(TR.y[int(i)])
        radius.append(min)
        radius2.append(min2)

#  print TR.sens
#  print peak
#  print amp
#  print radius
#  print radius2

  if len(peak)>0:
#    print len(peak)

    if len(peak)==1:
      sig=[radius[0]]
      popt,pcov=curve_fit(oneGauss2,x,n,p0=[peak[0],amp[0],sig[0]])
      if popt[1]>0:
        TR.popt.append(popt[0])
        TR.popt.append(popt[1])
        TR.popt.append(popt[2])

    elif len(peak)==2:
      sig=[radius[0],radius[1]]
      popt,pcov=curve_fit(twoGauss2,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1]])

      for i in range(6):
        TR.popt.append(popt[i])
 
    else:
      sig=[radius[0],radius[1],radius[2]]
      popt,pcov=curve_fit(threeGauss2,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1],peak[2],amp[2],sig[2]])
#      pb.plot(x,nGauss2(x,popt),color=col3,linewidth=1)
#      pb.show()

      for i in range(9):
        TR.popt.append(popt[i])

#  print TR.popt
  pb.plot(x,nGauss2(x,TR.popt),color=col3,linewidth=1)
#  pb.show()

  i=re.search('[0-9]+',filename)
  pb.savefig('figures/'+i.group()+'.png')
  pb.clf()

  return TR.popt

#-#-#-#-#-#-# 
#-#-#-#-#-#-# 

def Analyze(filename,SHOW=False,SAVE=True):
  ## start analysis
  TR=TRACE(filename) # initialize the trace object
  TR.get_scanparams() # get bias and pos and dist
  TR.getItrace() # add the data to TR.I
  TR.getHist() # get histogram
#  histFit(TR)
  smooth(TR)
  peak_finder(TR)
  threeGaussFit(TR)
  TR.params=hidden_peak(TR,filename)
  return 
  TR.resetHist()
  TR.plotHist(SHOW=SHOW,SAVE=SAVE)
  return
  #t1.smooth(TR)
  #t1.peak_finder(TR)
  #t1.threeGaussFit(TR)
  
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
