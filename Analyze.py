#!/usr/bin/python

import sys
import numpy as np

import ttCorr as ttC

numbin=300
Imax=100
Fs=1e4

def histFit(TR):
  guessparams=ttC.get_guessparams(TR)
  TR.params=ttC.three.threeGaussFit(TR,mode='vocal',guess=guessparams)
  ttC.showHist(TR)
  
  ### Fitting the histogram with Gaussians: 
  success=False
  while success==False:
    var=raw_input("OK or adjust guess ('y', '1neg', '1pos', '2neg', '2pos', or '3' for single, double, or triple gauss fit) ?")
    if var=='y':
      success=True
      tot=TR.params[0]+TR.params[2]+TR.params[4]
      totalCounts=tot*TR.numbin/TR.Imax
      TR.Pneg=TR.params[0]/tot #negative steady state
      TR.Pneu=TR.params[2]/tot #neutral steady state
      TR.Ppos=TR.params[4]/tot #positive steady state
      TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=TR.params   
      TR.SIGneg,TR.SIGneu,TR.SIGpos=(ttC.sigofmu(TR.params[1]),ttC.sigofmu(TR.params[3]),ttC.sigofmu(TR.params[5]))
    elif var=='3':
      guessparams=ttC.three.adjustGuess(TR)
      TR.params=ttC.three.threeGaussFit(TR,mode='vocal',guess=guessparams)
      ttC.showHist(TR)
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y':
          succ=True
          success=True
          TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=TR.params
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(ttC.sigofmu(TR.params[1]),ttC.sigofmu(TR.params[3]),ttC.sigofmu(TR.params[5]))
          tot=TR.params[0]+TR.params[2]+TR.params[4]
          totalCounts=tot*numbin/Imax
          TR.Pneg=TR.params[0]/tot #negative steady state
          TR.Pneu=TR.params[2]/tot #neutral steady state
          TR.Ppos=TR.params[4]/tot #positive steady state
        elif var=='n':
          succ=True
          success=False
        elif var!='':
          print "Ain't nobody got time for that!"
    elif var=='2neg':
      TR.histtype='2neg'
      params=params[:4]
      guessparams=ttC.two.adjustGuess(TR.x,TR.n,params)
      params=ttC.twoGaussFit(TR.x,TR.n,mode='vocal',guess=guessparams)
      ttC.showHist(TR.x,TR.n,tuple(params)+(0,-1000))
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y': 
          succ=True
          success=True
          TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=tuple(params)+(0.0,np.nan)
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(ttC.sigofmu(params[1]),ttC.sigofmu(params[3]),np.nan)
          tot=params[0]+params[2]
          totalCounts=tot*numbin/Imax
          TR.Pneg=params[0]/tot #negative steady state
          TR.Pneu=params[2]/tot #neutral steady state
          TR.Ppos=0
        elif var=='n':
          succ=True
          success=False
        elif var!='':
          print "Ain't nobody got time for that!"
    elif var=='2pos':
      TR.histtype='2pos'
      params=params[2:]
      guessparams=ttC.two.adjustGuess(TR.x,TR.n,params)
      params=ttC.twoGaussFit(TR.x,TR.n,mode='vocal',guess=guessparams)
      ttC.showHist(TR.x,TR.n,(0,-1000)+tuple(params))
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y':
          succ=True 
          success=True
          TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=(0.0,np.nan)+tuple(params)
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(np.nan,ttC.sigofmu(params[1]),ttC.sigofmu(params[3]))
          tot=params[0]+params[2]
          totalCounts=tot*numbin/Imax
          TR.Pneg=0 #negative steady state
          TR.Pneu=params[0]/tot #neutral steady state
          TR.Ppos=params[2]/tot
        elif var=='n':
          succ=True
          success=False
        elif var!='':
          print "Ain't nobody got time for that!"
    elif var=='1neg':
      TR.histtype='1neg'
      params=params[:2]
      guessparams=ttC.one.adjustGuess(TR.x,TR.n,params)
      params=ttC.one.oneGaussFit(TR.x,TR.n,mode='vocal',guess=guessparams)
      ttC.showHist(TR.x,TR.n,tuple(params)+(0,-1000,0,-1000))
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y': 
          succ=True
          success=True
          tot=params[0]
          totalCounts=tot*numbin/Imax
          TR.Pneg=1.
          TR.Pneu=0.
          TR.Ppos=0.
          TR.Aneg,TR.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=tuple(params)+(0.0,np.nan,0.0,np.nan)
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(ttC.sigofmu(params[1]),np.nan,np.nan)
        elif var==n:
          succ=True
          success=False
        elif var!='' and var!='n':
          print "Ain't nobody got time for that!"
    elif var=='1pos':
      states='1pos'
      params=params[4:]
      guessparams=ttC.one.adjustGuess(TR.x,TR.n,params)
      params=ttC.one.oneGaussFit(x,n,mode='vocal',guess=guessparams)
      ttC.showHist(TR.x,TR.n,(0,-1000,0,-1000)+tuple(params))
      succ=False
      while succ==False:
        var=raw_input("Happy yet ('y' or 'n') ?")
        if var=='y': 
          succ=True
          success=True
          tot=params[0]
          totalCounts=tot*numbin/Imax
          TR.Pneg=0.
          TR.Pneu=0.
          TR.Ppos=1.
          TR.Aneg,Tr.MUneg,TR.Aneu,TR.MUneu,TR.Apos,TR.MUpos=(0.0,np.nan,0.0,np.nan)+tuple(params)
          TR.SIGneg,TR.SIGneu,TR.SIGpos=(np.nan,np.nan,ttC.sigofmu(params[1]))
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

######
######
######
######

def Analyze(filename):
  TR=ttC.TRACE(filename) # initialize the trace object
  TR=ttC.get_scanparams(TR) # get bias and pos and dist
  TR=ttC.getItrace(TR) # add the data to TR.I
  TR=ttC.getHist(TR) # get histogram
  TR=histFit(TR) # fit the histogram
  
  taumax=10e-3 # 10 ms
  if TR.histtype=='3':
    # get time correlations over three ranges
    IrangeNEG,IrangeNEU,IrangePOS,PoALL=ttC.three.ranges(TR.params)
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
    NEG=ttC.three.getCorr(TR,IrangeNEG,taumax) # NEG=(T,P1,P2,P3)
    print '\nWorking on neutral range:'
    NEU=ttC.three.getCorr(TR,IrangeNEU,taumax) # NEU=(T,P1,P2,P3)
    print '\nWorking on positive range:'
    POS=ttC.three.getCorr(TR,IrangePOS,taumax) # POS=(T,P1,P2,P3)
    print '\n'
  elif TR.histtype=='2neg':
    # get time correlations over three ranges
    IrangeONE,IrangeTWO,PoALL=ttC.two.ranges(TR.params)
    gamma=TR.Pneg/TR.Pneu #this gives f from e... f=gamma*e
    q=(gamma,)
    print '\nWorking on first range:'
    NEG=ttC.two.getCorr(TR,IrangeONE,taumax) # NEG=(T,P1,P2,P3)
    print '\nWorking on second range:'
    NEU=ttC.two.getCorr(TR,IrangeTWO,taumax) # NEU=(T,P1,P2,P3)
    #print '\nWorking on positive range:'
    #POS=ttC.getCorr2(I,IrangePOS,taumax,params) # POS=(T,P1,P2,P3)
    print '\n'
  elif TR.histype=='2pos':
    #print '\nWorking on negative range:'
    #NEG=ttC.getCorr2(I,IrangeNEG,taumax,params) # NEG=(T,P1,P2,P3)
    IrangeONE,IrangeTWO,PoALL=ttC.two.ranges(TR.params)
    gamma=TR.Pneu/TR.Ppos #this gives f from e... f=gamma*e
    q=(gamma,)
    print '\nWorking on neutral range:'
    NEU=ttC.two.getCorr(TR,IrangeONE,taumax) # NEU=(T,P1,P2,P3)
    print '\nWorking on positive range:'
    POS=ttC.two.getCorr(TR,IrangeTWO,taumax) # POS=(T,P1,P2,P3)
    print '\n'
  elif TR.histtype=='1neg' or TR.histtype=='1pos':
    print "JUST ONE STATE, YO!"
  
  #fitting:
  if TR.histtype=='3':
    print q
    print PoALL
    p=ttC.three.dublexFit(np.array(NEG[0]),np.array(NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]),q,PoALL)
    a,b=p
    TR.a=a#*TR.Fs # converts to units of Hz --- no longer necessary
    TR.b=b#*TR.Fs # converts to units of Hz --- no longer necessary
    TR.c,TR.d=(a*alpha,b*beta)
    print 'a,b,c,d are:'
    print str('%.4f, '%TR.a)+str('%.4f, '%TR.b)+str('%.4f, '%TR.c)+str('%.4f, '%TR.d)
    ttC.three.dublexPLOT(NEG,NEU,POS,p,q,PoALL,taumax)
  elif TR.histtype=='2neg':
    p=ttC.two.singlexFit(np.array(NEG[0]),np.array(NEG[1]+NEG[2]+NEU[1]+NEU[2]),q,PoALL)
    a,=p
    TR.a=a*TR.Fs # converts to units of Hz
    TR.c=TR.a*gamma
    TR.b=np.nan
    TR.d=np.nan
    print 'a,b,c,d are:'
    print str('%.4f, '%TR.a)+str('%.4f, '%TR.b)+str('%.4f, '%TR.c)+str('%.4f, '%TR.d)
    #ttC.two.singlexPLOT(NEG,NEU,p,q,PoALL,taumax)
    ttC.two.singlexPLOT(NEG,NEU,p,q,PoALL,taumax)
  elif TR.histtype=='2pos':
    p=ttC.two.singlexFit(np.array(NEU[0]),np.array(NEU[1]+NEU[2]+POS[1]+POS[2]),q,PoALL)
    b,=p
    TR.b=b*TR.Fs # converts to units of Hz
    TR.d=TR.b*gamma
    TR.a=np.nan
    TR.c=np.nan
    print 'a,b,c,d are:'
    print str('%.4f, '%TR.a)+str('%.4f, '%TR.b)+str('%.4f, '%TR.c)+str('%.4f, '%TR.d)
    #ttC.two.singlexPLOT(NEU,POS,p,q,PoALL,taumax)
    ttC.two.singlexPLOT(NEU,POS,p,q,PoALL,taumax)
  elif TR.histtype=='1neg' or TR.histtype=='1pos':
    TR.a=np.nan
    TR.b=np.nan
    TR.c=np.nan
    TR.d=np.nan
    print 'a,b,c,d are:'
    print str('%.4f, '%TR.a)+str('%.4f, '%TR.b)+str('%.4f, '%TR.c)+str('%.4f, '%TR.d)
  
  return TR

####################
# Main Function. 
####################

def main():
  if len(sys.argv)==1: # no arguments passed by user. must stop
    print '\n/-----\\\n|ERROR|===>>> Try inputs of the form python Analyze.py foo.dat (mode)\n\-----/\n'
    sys.exit(0) # leave the program
  filename=sys.argv[1]
  if len(sys.argv)==3:  
    MODE=sys.argv[2] # modes are meant to produce a particular plot or particular data... ALL gets them all
  else:
    MODE=None
  if MODE in ['histogram','Irange','trace','traceCorr','corrHists',None]:
    print '\nRunning Analyze.py in mode "'+str(MODE)+'"\n'
  else:
    print "\nWatch it son!\nYou are only allowed modes in:\n['histogram','Irange','trace',traceCorr,corrHists,None]\n"
    sys.exit(0)
  Analyze(filename)
  return

if __name__=='__main__': # trick to allow using main as a function
  main()                 # main function is called here
