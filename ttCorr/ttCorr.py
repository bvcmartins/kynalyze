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
#import threestate as three
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
  def __init__(self,fname,numbin=300,Imax=None,Fs=1e4,taumax=10e-3):
    
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
    self.ranges=None
    self.taumax=taumax
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
######
    ## based on subsequent methods
    self.get_scanparams()
    self.getItrace()
    if self.Imax==None: self.Imax=max(self.I)

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
      self.I=[float(i)*1E12 for i in I]    #change units to pA
    else: #assume that it's one of my simple data files
      data=f.readlines()
      self.I=[float(dt.split()[1])*1E12 for dt in data]
    return None

  #make histogram data
  def getHist(self):
    n,bins=np.histogram(self.I,bins=np.linspace(0,self.Imax,self.numbin+1))
    self.n=n
    self._n=copy.deepcopy(n) # this is a backup so that we can reset n later. It means we can smooth or alter n as much as we want and still return to the raw histogram.
    self.bins=bins
    self.x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
    return
  
  #fit histogram
  def fitHist(self):
    smooth(self)
    peak_finder(self)
    threeGaussFit(self)
    self.params=hidden_peak(self)
    self.resetHist()
    return
  
  #cut out bad peaks (replace this with something more sophisticated later)
  def peakfilter(self):
    i=0
    ## cut out bad or insignificant peaks
    while i < len(self.params)/3:
      if self.params[i*3+1]<0: self.params[i*3:i*3+3]=[]
      elif self.params[i*3+2]<0: self.params[i*3:i*3+3]=[]
      elif self.params[i*3+1]/(len(self.I)*self.Imax/self.numbin)<0.01: self.params[i*3:i*3+3]=[]
      else: i+=1
    ## order peaks
    temp=[]
    for i in range(len(self.params)/3): temp.append(self.params[3*i:3*i+3])
    temp.sort(key=lambda el:el[0])
    self.params=[]
    for i in range(len(temp)): self.params+=temp[i]
    return
  
  def getRanges(self):
    self.ranges=[]
    for i in range(len(self.params)/3):
      self.ranges.append((self.params[3*i]-self.params[3*i+2],self.params[3*i]+self.params[3*i+2])) ## (Imin,Imax)
    return
  
  #plot histogram
  def plotHist(self,pars=None,SHOW=True,SAVE=True): 
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
        pb.plot(np.linspace(0,100,500),nGauss(np.linspace(0,100,500),*params[i*3:i*3+3]),color=colours[i],linewidth=1)
      pb.plot(np.linspace(0,100,500),nGauss(np.linspace(0,100,500),*params),'--k',linewidth=1)
    if self.ranges!=None:
      for i in range(len(self.ranges)):
        rn=self.ranges[i]
        mask=(self.x>rn[0])&(self.x<rn[1])
        pb.bar(np.array(self.x)[mask],np.array(self.n)[mask],float(self.Imax)/self.numbin,color=colours[i],linewidth=0.4,align='center')
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
    else: self.plotHist(SHOW=SHOW) # this looks problematic
    #if savename!=None:pb.savefig(savename)
    if SHOW:pb.show()
    return
    
  def getCorr(self):
    Tmax=int(self.taumax*self.Fs) ## integer corresponding to the index
    
    #-#-# define the constrained gaussian function:
    def nGauss_constrained(x,*amplitudes):
      G=[]
      for i in range(len(amplitudes)):
        mu,A,sig=(self.params[3*i],amplitudes[i],self.params[3*i+2])
        G.append((A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig)))
      return sum(G)
    #-#-#
    
    self.Pcorr=np.zeros([len(self.ranges),len(self.params)/3,Tmax+1]) # Pcorr=[rangenumber,statenumber,T]
    self.PcorrERR=np.zeros([len(self.ranges),len(self.params)/3,Tmax+1])
    #self.ERRcorr=np.zeros([len(self.ranges),Tmax+1]) # Pcorr=[rangenumber,statenumber,T]
    for irn in range(len(self.ranges)):
      rn=self.ranges[irn]
      for T in range(Tmax+1):
        indeces=np.arange(len(self.I[:(-T)]))[(self.I[:(-T)]>rn[0])&(self.I[:(-T)]<rn[1])]+T
        Icorr=np.array(self.I)[indeces]
        nCorr,binsCorr=np.histogram(Icorr,bins=self.bins)
        xCorr=[0.5*(binsCorr[k]+binsCorr[k+1]) for k in range(len(binsCorr)-1)]
        amplitudes=self.params[1::3]
        popt,pcov=curve_fit(nGauss_constrained,xCorr,nCorr,p0=amplitudes)
        if not ((np.max(pcov)==np.inf)|(np.min(pcov)==-np.inf)): ## this ignores the first ugly data point, where the fits are bad
          self.Pcorr[irn,:,T]=np.array(popt)/popt.sum()
          self.PcorrERR[irn,:,T]=np.sqrt(np.diagonal(pcov))/popt.sum()
        else:
          self.Pcorr[irn,:,T]=np.nan
          self.PcorrERR[irn,:,T]=np.nan
    return
    
  def resetHist(self):
    self.n=copy.deepcopy(self._n)
    return
  
  #plot time-time correlations
  def plotCorr(self,pars=None,SHOW=True,SAVE=True): 
    skipfits=True
    pb.clf()
    for irn in range(len(self.ranges)):
      for ist in range(len(self.params)/3):
        pb.errorbar(np.arange(len(self.Pcorr[irn,ist,:]))/self.Fs,self.Pcorr[irn,ist,:],yerr=self.PcorrERR[irn,ist,:],color=colours[ist])
      pb.ylim([0,1])
      pb.grid(True)
      pb.xlabel('Time (s)',fontsize=20)
      pb.ylabel('Probabilities',fontsize=20)
      pb.xticks(fontsize=16)
      pb.yticks(fontsize=16)
      if SAVE:
        fn=os.path.basename(self.filename)
        pb.savefig(self.figdir+'corr_'+'range_'+str(irn)+fn[:-4]+'.eps')
        f=open(self.datdir+fn[:-4]+'.dat','a')
        f.seek(0,2)
        f.write('#corr_range_'+str(irn)+'\n')
        #f.write('BinCentre(pA)\tBinMin(pA)\tBinMax(pA)\tCounts\n')
        #for i in range(len(self.n)):
        #  f.write(str(self.x[i])+'\t'+str(self.bins[i])+'\t'+str(self.bins[i+1])+'\t'+str(self.n[i])+'\n')
        f.close()
      if SHOW:pb.show()
      pb.clf()
    return
  
####

#-#-# This needs to be replaced with something much more flexible:
#-#-#

def sigofmu(mu):
  sig=0.5+(1./20)*mu
  return sig

#-#-# flexible gaussian fit:
#-#-#

def nGauss(x,*params):
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

def peak_finder(TR):
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]
  for i in np.where(pk==1)[0]:
    TR.peak.append(TR.x[int(i)])
    TR.amp.append(TR.y[int(i)])
  pb.bar(TR.x,TR.n,float(TR.Imax)/TR.numbin,color='0.8',linewidth=0.4,align='center')

def threeGaussFit(TR):

  if len(TR.peak)>0:
    if len(TR.peak)==1:
      sig=[sigofmu(TR.peak[0])]
      popt,pcov=curve_fit(nGauss,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0]])
      TR.popt.append(popt[0])
      TR.popt.append(popt[1])
      TR.popt.append(popt[2])
      TR.n=TR.n-nGauss(TR.x,popt[0],popt[1],popt[2])

    elif len(TR.peak)==2:
      sig=[sigofmu(TR.peak[0]),sigofmu(TR.peak[1])]
      popt,pcov=curve_fit(nGauss,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1]])

      for i in range(6):
        TR.popt.append(popt[i])
      TR.n=TR.n-nGauss(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

    else:
      sig=[sigofmu(TR.peak[0]),sigofmu(TR.peak[1]),sigofmu(TR.peak[2])]
      popt,pcov=curve_fit(nGauss,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1],TR.peak[2],TR.amp[2],sig[2]])
      for i in range(9):
        TR.popt.append(popt[i])
      TR.n=TR.n-nGauss(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8])

def hidden_peak(TR):
  window_len=11
  window='hanning'
  s=np.r_[TR.n[window_len-1:0:-1],TR.n,TR.n[-1:-window_len:-1]]
  w=eval('np.'+window+'(window_len)')
  TR.y=np.convolve(w/w.sum(),s,mode='valid')
  TR.y=TR.y[5:-5]
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]

  x=TR.x
  n=TR.n
 
  peak=[]
  amp=[]

  for i in np.where(pk==1)[0]:
    if TR.y[int(i)]>TR.sens:
      peak.append(TR.x[int(i)])
      amp.append(TR.y[int(i)])

  if len(peak)>0:

    if len(peak)==1:
      sig=[sigofmu(peak[0])]
      popt,pcov=curve_fit(nGauss,x,n,p0=[peak[0],amp[0],sig[0]])
      TR.popt.append(popt[0])
      TR.popt.append(popt[1])
      TR.popt.append(popt[2])

    elif len(peak)==2:
      sig=[sigofmu(peak[0]),sigofmu(peak[1])]
      popt,pcov=curve_fit(nGauss,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1]])

      for i in range(6):
        TR.popt.append(popt[i])
    
    else:
      sig=[sigofmu(peak[0]),sigofmu(peak[1]),sigofmu(peak[2])]
      popt,pcov=curve_fit(nGauss,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1],peak[2],amp[2],sig[2]])
 
      for i in range(9):
        TR.popt.append(popt[i])
  #i=re.search('[0-9]+',TR.filename)
  return TR.popt

#-#-#-#-#-#-# 
#-#-#-#-#-#-# 

def Analyze(filename,SHOW=False,SAVE=True):
  ## start analysis
  TR=TRACE(filename,numbin=300,Imax=100) # initialize the trace object, get parameters, and read in the trace
  TR.getHist() ## get histogram - set TR.bins, TR.x, TR.n
  TR.fitHist() ## fit histogram - set TR.params
  TR.peakfilter()
  TR.getRanges() ## set TR.ranges
  TR.plotHist(SHOW=SHOW,SAVE=SAVE) ## this is a slow step right now... but that may not be a big deal.
  TR.getCorr() ## sets TR.Pcorr - a 3D array of conditional probabilities: Pcorr[initialrange,to_state,TAU]
  TR.plotCorr(SHOW=SHOW,SAVE=SAVE) ## makes plots from Pcorr
  return
  
  ## get time correlations over three ranges
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
