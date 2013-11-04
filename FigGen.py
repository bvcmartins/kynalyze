import os
import scipy.fftpack
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import math

import sys
import numpy as np
import pylab as pb

import ttCorr as ttC
import Analyze
from matplotlib.gridspec import GridSpec
import pickle

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

numbin=300
Imax=100
Fs=1e4

def MTprint(MESSAGE):
  print ' '*len(MESSAGE)
  print '-'*len(MESSAGE)
  print MESSAGE
  print '-'*len(MESSAGE)
  print ' '*len(MESSAGE)
  return

def Exp(p,x):
  A,lam = p
  return A*np.exp(np.array(x)/lam)

def ExpErr(p, x, y): 
  error=Exp(p,x) - y
  return error

def ExpFit(x,y,guess=(1.,1.)):
  params,success=leastsq(ExpErr,guess,args=(x,y))
  A,lam=params
  return params

def Analyze(TR,TAU):
  
  # get and fit histogram of data
  states='3'
  guessparams=ttC.get_guessparams(TR)
  TR.params=ttC.three.threeGaussFit(TR,mode='vocal',guess=guessparams)
  
  ### Fitting the histogram with Gaussians:
  MTprint("I'm going to assume that the triple Gaussian fit is good right off the bat.")
  tot=TR.params[0]+TR.params[2]+TR.params[4]
  totalCounts=tot*numbin/Imax
  Pneg=TR.params[0]/tot #negative steady state
  Pneu=TR.params[2]/tot #neutral steady state
  Ppos=TR.params[4]/tot #positive steady state
  print 'Steady state probabilities (from multiple Gaussian fit):'
  print 'Pneg: '+str('%.4f'%Pneg)
  print 'Pneu: '+str('%.4f'%Pneu)
  print 'Ppos: '+str('%.4f'%Ppos)+'\n'
  if np.abs(totalCounts-len(TR.I))>(float(len(TR.I))/100):
    MTprint('---> WARNING: totalCounts is ' +str(totalCounts))

  taumax=10.0
  # get time correlations over three ranges
  IrangeNEG,IrangeNEU,IrangePOS,PoALL=ttC.three.ranges(TR.params)
  alpha=Pneg/Pneu
  beta=Pneu/Ppos
  q=(alpha,beta)
  print '\nWorking on negative range:'
  NEG=ttC.three.getCorr(TR,IrangeNEG,taumax) # NEG=(T,P1,P2,P3)
  print '\nWorking on neutral range:'
  NEU=ttC.three.getCorr(TR,IrangeNEU,taumax) # NEU=(T,P1,P2,P3)
  print '\nWorking on positive range:'
  POS=ttC.three.getCorr(TR,IrangePOS,taumax) # POS=(T,P1,P2,P3)
  print '\n'
  
  
  p=ttC.three.dublexFit(np.array(NEG[0]),np.array(NEG[1]+NEG[2]+NEG[3]+NEU[1]+NEU[2]+NEU[3]+POS[1]+POS[2]+POS[3]),q,PoALL)
  a,b=p
  a=a*Fs # converts to units of Hz
  b=b*Fs # converts to units of Hz
  c,d=(a*alpha,b*beta)
  print 'a,b,c,d are:'
  print str('%.4f, '%a)+str('%.4f, '%b)+str('%.4f, '%c)+str('%.4f, '%d)
  dublexPLOT(NEG,NEU,POS,p,q,PoALL,taumax,TAU)
  
  return

def dublexPLOT(NEG,NEU,POS,p,q,PoALL,taumax,TAU):
#NEG NEU POS p q taumax
  T=NEG[0]
  tau=[float(t)*1000/Fs for t in T] #convert to ms
  TMAX=taumax*1000/Fs #convert to ms
  
  subsize=0.39
  C1=0.09
  C2=0.58
  R1=0.57
  R2=0.08
  
  dataopts={'linestyle':'None','linewidth':0.3,'marker':'.','markersize':1}
  gs2=pb.GridSpec(1,1)
  gs2.update(left=C2, right=C2+subsize, bottom=R1,top=R1+subsize)
  ax=pb.subplot(gs2[0,0])
  #ax2=pb.subplot2grid((4,4),(0,2),rowspan=2,colspan=2)
  pNEG1,=pb.plot(tau,NEG[1],col1,**dataopts)
  pNEG2,=pb.plot(tau,NEG[2],col2,**dataopts)
  pNEG3,=pb.plot(tau,NEG[3],col3,**dataopts)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[0][1],col1)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[0][2],col2)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[0][3],col3)
  pb.text(0.5,0.9,'(e)',fontsize=10)
  for tau_ in TAU[1:]:
    pb.plot([tau_,tau_],[0,1],'--k',linewidth=1)
    pb.text(tau_+0.1,0.5,r"$\tau$"+" = %.1f"%tau_+"ms",fontsize=6,rotation=-90)
  lg=pb.legend([pNEG1, pNEG2, pNEG3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=6,markerscale=2,borderaxespad=0.9)
  fr = lg.get_frame()
  fr.set_lw(0.2)
  [line.set_zorder(3) for line in ax.lines]
  pb.grid(True,color='0.7')
  pb.xlim([0,TMAX])
  pb.ylim([0,1])
  pb.xticks(fontsize=6)
  pb.yticks(fontsize=6)
  pb.xlabel(r"$\tau$ (ms)",fontsize=8,labelpad=1)
  pb.ylabel('P(t)',fontsize=8,labelpad=1)
  
  gs3=pb.GridSpec(1,1)
  gs3.update(left=C1, right=C1+subsize, bottom=R2,top=R2+subsize)
  ax=pb.subplot(gs3[0,0])
  #pb.subplot2grid((4,4),(2,0),rowspan=2,colspan=2)
  pb.grid(True,color='0.7')
  pNEU1,=pb.plot(tau,NEU[1],col1,**dataopts)
  pNEU2,=pb.plot(tau,NEU[2],col2,**dataopts)
  pNEU3,=pb.plot(tau,NEU[3],col3,**dataopts)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[1][1],col1)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[1][2],col2)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[1][3],col3)
  pb.text(0.5,0.9,'(f)',fontsize=10)
  lg=pb.legend([pNEU1, pNEU2, pNEU3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=6,markerscale=2,borderaxespad=0.9)
  fr = lg.get_frame()
  fr.set_lw(0.2)
  [line.set_zorder(3) for line in ax.lines]
  pb.xlim([0,TMAX])
  pb.ylim([0,1])
  pb.xticks(fontsize=6)
  pb.yticks(fontsize=6)
  pb.xlabel(r"$\tau$ (ms)",fontsize=8,labelpad=1)
  pb.ylabel('P(t)',fontsize=8,labelpad=1)
  
  gs4=pb.GridSpec(1,1)
  gs4.update(left=C2, right=C2+subsize, bottom=R2,top=R2+subsize)
  ax=pb.subplot(gs4[0,0])
  #pb.subplot2grid((4,4),(2,2),rowspan=2,colspan=2)
  pPOS1,=pb.plot(tau,POS[1],col1,**dataopts)
  pPOS2,=pb.plot(tau,POS[2],col2,**dataopts)
  pPOS3,=pb.plot(tau,POS[3],col3,**dataopts)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[2][1],col1)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[2][2],col2)
  pb.plot(tau,ttC.three.dublex(p,T,q,PoALL)[2][3],col3)
  pb.text(0.5,0.9,'(g)',fontsize=10)
  lg=pb.legend([pPOS1, pPOS2, pPOS3], ["$P_{-}$", "$P_{o}$", "$P_{+}$"],fontsize=6,markerscale=2,borderaxespad=0.9)
  fr = lg.get_frame()
  fr.set_lw(0.2)
  [line.set_zorder(3) for line in ax.lines]
  pb.grid(True,color='0.7')
  pb.xlim([0,TMAX])
  pb.ylim([0,1])
  pb.xticks(fontsize=6)
  pb.yticks(fontsize=6)
  pb.xlabel(r"$\tau$ (ms)",fontsize=8,labelpad=1)
  pb.ylabel('P(t)',fontsize=8,labelpad=1)
  
  #if directory !=None: pb.savefig(directory+'POS.png')
  #if directory!=None: return
  #else:
  #  pb.show()
  return

def corrHist(I,Irange,TAU):
  subsize=0.39
  C1=0.09
  C2=0.58
  R1=0.57
  R2=0.08
  gs1=pb.GridSpec(2,2)
  gs1.update(left=C1, right=C1+subsize, bottom=R1,top=R1+subsize, wspace=0.06, hspace=0.06)
  letter=['a','b','c','d']
  i=1
  for tau in TAU:
    print (int((i-1)/2),(i-1)%2)
    ax=pb.subplot(gs1[int((i-1)/2),(i-1)%2])
    n,bins=ttC.getHist(I)
    pb.bar(bins[1:],n,float(Imax)/numbin,color='0.5',alpha=0.5,linewidth=0)
    nCorr,binsCorr=ttC.corrHist(I,Irange,int(tau*Fs/1000))
    pb.bar(binsCorr[1:],nCorr,float(Imax)/numbin,color=col1,linewidth=0)
    pb.xlim([10,35])
    pb.ylim([0,1.2*max(n)])
    #pb.grid(True)
    #pb.xlabel('I (pA)',fontsize=16)
    #pb.ylabel('counts',fontsize=16)
    #pb.xticks([])
    #pb.yticks([])
    pb.text(10.3,0.98*max(n),'('+letter[i-1]+')',fontsize=8)
    pb.text(34,1.16*max(n),(r"$\tau$"+" = %.1f"%tau)+'ms',fontsize=6,ha='right',va='top')
    if i==1:
      pb.ylabel('counts',fontsize=8)
      pb.xticks([])
      pb.yticks([])
    elif i==2:
      pb.xticks([])
      pb.yticks([])
    elif i==3:
      pb.xlabel('$I_T$ (pA)',fontsize=8,labelpad=1)
      pb.ylabel('counts',fontsize=8)
      pb.yticks([])
      pb.xticks([10,20,30],fontsize=6)
      ax.xaxis.set_ticks_position('bottom')
    elif i==4:
      pb.xlabel('$I_T$ (pA)',fontsize=8,labelpad=1)
      pb.yticks([])
      pb.xticks([10,20,30],fontsize=6)
      ax.xaxis.set_ticks_position('bottom')
    i+=1
  return
  
def showTrace(TR):
  t=[float(i)*1000/TR.Fs for i in range(len(TR.I))]
  pb.clf()
  pb.plot(t,TR.I,'k',marker='.',markersize=0.2,linewidth=0.1)
  pb.grid(True)
  pb.xlim([200,400])
  pb.ylim([10,35])
  pb.xlabel('$t$(ms)',fontsize=8)
  pb.ylabel('$I_T$(pA)',fontsize=8)
  pb.xticks(fontsize=6)
  pb.yticks(fontsize=6)
  pb.tight_layout(pad=0.1,rect=(0,0,0.99,0.97))
  return

#show histogram
def showHist(x,n,params,xlim=[0,60]):
  pb.bar(x,n,float(Imax)/numbin,color='0.9',linewidth=0.05,align='center')
  pb.plot(np.linspace(0,100,500),ttC.oneGauss(params[0:2],np.linspace(0,100,500)),color=col1,linewidth=0.6)
  pb.plot(np.linspace(0,100,500),ttC.oneGauss(params[2:4],np.linspace(0,100,500)),color=col2,linewidth=0.6)
  pb.plot(np.linspace(0,100,500),ttC.oneGauss(params[4:6],np.linspace(0,100,500)),color=col3,linewidth=0.6)
  pb.plot(np.linspace(0,100,500),ttC.threeGauss(params,np.linspace(0,100,500)),'--k',linewidth=0.8)
  pb.xlim([10,35])
  pb.grid(True)
  pb.xlabel('$I_T$(pA)',fontsize=8)
  pb.ylabel('Counts',fontsize=8)
  pb.xticks(fontsize=6)
  pb.yticks(fontsize=6)
  pb.tight_layout(pad=0.1,rect=(0,0,0.99,0.97))
  return

def cmap(volt,letter):
  basename='./data/telegraph-const-height-'+volt+'-'
  x=pb.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
  x=[i*10**9 for i in x]
  y=pb.arange(numbin)*Imax/numbin
  X,Y=pb.meshgrid(x,y)
  H=pb.zeros(pb.shape(X))
  for num in range(25)[1:]:
    print num
    filename=basename+('%(number)03d' % {'number': num})+'.dat'
    TR=ttC.TRACE(filename,Fs=1e4)
    TR=ttC.getItrace(TR)
    #tup=pb.hist(I,bins=numbin,range=(0,100),color='b')
    n,bins=np.histogram(TR.I,bins=np.linspace(0,TR.Imax,TR.numbin+1))
    H[:,num-1]=pb.array(n)
  #pb.clf()
  pb.pcolor(X,Y,H,cmap='hot')
  #pb.colorbar()
  pb.xlim([min(x),max(x)])
  pb.ylim([0,99])
  if letter=='e':
    pb.xticks([2,3,4],fontsize=14)
    pb.yticks(fontsize=14)
    pb.xlabel('d (nm)',fontsize=18)
    pb.ylabel('I (pA)',fontsize=18)
  elif letter=='a':
    pb.yticks(fontsize=14)
    pb.ylabel('I (pA)',fontsize=18)
    pb.xticks([])
  elif letter=='f' or letter=='g' or letter=='h':
    pb.xticks([2,3,4],fontsize=14)
    pb.xlabel('d (nm)',fontsize=18)
    pb.yticks([])
  else:
    pb.xticks([])
    pb.yticks([])
  pb.text(3.75,90,volt[:1]+'.'+volt[-2:]+'V',fontsize=14,color='w')
  pb.text(2.1,87,'('+letter+')',fontsize=18,color='w')
  #pb.savefig('./cmaps/image-'+volt+'.png')
  #pb.show()
  #pb.close()
  return
  
def cmaps():
  volts=['1pt'+num for num in ['30','35','40','45','50','55','60','65']]
  letters=['a','b','c','d','e','f','g','h']
  #pb.figure(1,figsize=(12,6))
  pb.figure(1,figsize=(12,5))
  pb.clf()
  i=1
  for volt in volts:
    print volt+" votls"
    pb.subplot(2,4,i)
    cmap(volt,letters[i-1])
    i+=1
  pb.tight_layout()
  pb.savefig('./cmaps.png',dpi=900)
  pb.show()
  return
  
def onecmap(basename='./data/telegraph-const-height-1pt45-'):
  #basename='./data/telegraph-const-height-1pt45-'
  x=pb.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
  x=[i*10**9 for i in x]
  y=pb.arange(numbin)*Imax/numbin
  X,Y=pb.meshgrid(x,y)
  H=pb.zeros(pb.shape(X))
  for num in range(25)[1:]:
    filename=basename+('%(number)03d' % {'number': num})+'.dat'
    TR=ttC.TRACE(filename,Fs=1e4)
    TR=ttC.getItrace(TR)
    #tup=pb.hist(I,bins=numbin,range=(0,100),color='b')
    n,bins=np.histogram(TR.I,bins=np.linspace(0,TR.Imax,TR.numbin+1))
    H[:,num-1]=pb.array(n)
  pb.pcolor(X,Y,H,cmap='hot')
  cbar=pb.colorbar(ticks=[0,800,1600,2400])
  cbar.set_label('counts',fontsize=8)
  cbar.ax.tick_params(labelsize=6) 
  pb.xlim([min(x),max(x)])
  pb.ylim([0,60])
  pb.xticks([2,3,4],fontsize=6)
  pb.yticks(fontsize=6)
  pb.xlabel('d (nm)',fontsize=8)
  pb.ylabel('I (pA)',fontsize=8)
  return
  
def main():
  if len(sys.argv)==1:
    print '\n/-----\\\n|ERROR|===>>> Try inputs of the form python Analyze.py MODE (optionally more) \n\-----/\n'
    sys.exit(0)
  MODE=sys.argv[1]
  if MODE not in ['cmaps','onecmap','threecmaps','histogram','trace','corr','rates']:
    MTprint("That is not a valid mode! Try entering cmaps, onecmap, threesmaps, histogram, trace, corr, or rates as modes")
    #print ' '*len(MESSAGE);print '-'*len(MESSAGE);print MESSAGE;print '-'*len(MESSAGE);print ' '*len(MESSAGE)
    sys.exit(0)
  
  if MODE=='cmaps':
    cmaps()
  
  elif MODE=='onecmap':
    pb.figure(1,figsize=(3.365,3.365*2/3))
    pb.clf()
    onecmap()
    pb.tight_layout()
    #cmap(volt,letters[i-1])
    pb.savefig('./onecmap.png',dpi=900,bbox_inches='tight')
    pb.show()
  
  elif MODE=='threecmaps':
    basenames=['./data/telegraph-const-height-1pt40-','./data/telegraph-const-height-1pt45-','./data/telegraph-const-height-1pt50-']
    letter=['a','b','c']
    pb.figure(1,figsize=(3.365,3.365*1.5))
    pb.clf()
    i=1
    for basename in basenames:
      pb.subplot(3,1,i)
      onecmap(basename)
      pb.text(3.9,52,basename[-6]+'.'+basename[-3:-1]+'V',fontsize=8,color='w')
      pb.text(2.1,50,'('+letter[i-1]+')',fontsize=12,color='w')
      i+=1
    pb.tight_layout()
    pb.savefig('./threecmaps.png',dpi=900,bbox_inches='tight')
    pb.show()
  
  elif MODE=='histogram':
    if len(sys.argv)!=3:
      MTprint("histogram mode requires you to specify a filename")
      sys.exit(0)
    filename=sys.argv[2]
    I=ttC.getItrace(filename)
    n,bins=ttC.getHist(I)
    x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
    params=ttC.threeGaussFit(x,n)
    pb.figure(1,figsize=(3.375/2,3.375/2))#,dpi=600)
    pb.clf()
    showHist(x,n,params)
    pb.savefig('histogram.eps')
    pb.show()
  
  elif MODE=='trace':
    if len(sys.argv)!=3:
      MTprint("histogram mode requires you to specify a filename")
      sys.exit(0)
    filename=sys.argv[2]
    TR=ttC.TRACE(filename)
    TR=ttC.getItrace(TR)
    pb.figure(1,figsize=(3.375,3.375/3))#,dpi=600)
    pb.clf()
    showTrace(TR)
    pb.savefig('trace.eps')
    pb.show()
  
  elif MODE=='corr':
    if len(sys.argv)!=3:
      MTprint("corr mode requires you to specify a filename")
      sys.exit(0)
    filename=sys.argv[2]
    TR=ttC.TRACE(filename)
    TR=ttC.getItrace(TR)
    TR=ttC.getHist(TR)
    #x=[0.5*(bins[k]+bins[k+1]) for k in range(len(bins)-1)] #center of bins
    TR.params=ttC.threeGaussFit(TR.x,TR.n)
    Irange=ttC.three.ranges(TR.params)[0]
    TAU=[0,0.5,2,8]#length should be 4
    pb.rcParams['xtick.major.pad']='3'
    pb.rcParams['ytick.major.pad']='3'
    pb.figure(1,figsize=(3.375,3.375))
    pb.clf()
    Analyze(TR,TAU)
    corrHist(TR.I,Irange,TAU)
    pb.savefig('./corr.eps')
    pb.show()
  
  elif MODE=='rates':
    
    COOL='b'#'g'#'#669966'
    WARM='g'#'#FFFF66'
    HOT='r'#'#FF0000'
    
    picklenames=['./PICKLES/1pt40_.p','./PICKLES/1pt45_.p','./PICKLES/1pt50_.p'] # aliases are A,B,C
    SET=pickle.load(open( './PICKLES/1pt40.p', 'r'))
    TUPLE=[{},{},{}]
    i=0
    for picklename in picklenames:
      print picklename
      SET=pickle.load(open( picklename, 'r'))
      dist=[]
      a=[];b=[];c=[];d=[]
      for TR in SET:
        if 1==1:#SET.V==1.40:
          dist.append(TR.dist)
          a.append(TR.a)
          b.append(TR.b)
          c.append(TR.c)
          d.append(TR.d)
      dict={'dist':dist,'a':a,'b':b,'c':c,'d':d}
      TUPLE[i]=dict
      i+=1
    A,B,C=TUPLE
    
    #FILLING RATES
    X=A['dist']+A['dist']+B['dist']+B['dist']+C['dist']+C['dist']
    Y=A['c']+A['d']+B['c']+B['d']+C['c']+C['d']
    poppers=[]
    for i in range(len(Y)):
      if math.isnan(Y[i]):
        poppers.append(i)
    poppers.reverse()
    for popper in poppers:
      X.pop(popper)
      Y.pop(popper)
    Amp,lam=ExpFit(X,Y,guess=(100000,-0.6))
    pb.figure(1,figsize=(3.375,3.375*5/3))
    pb.clf()
    pb.subplot(2,1,1)
    pb.plot(np.linspace(2,4.5,300),Exp((Amp,lam),np.linspace(2,4.5,300)),'--k')
    PLTAc,=pb.plot(A['dist'],A['c'],color=COOL,linestyle='None',marker='o',markersize=3,linewidth=0.1)
    PLTAd,=pb.plot(A['dist'],A['d'],color=COOL,linestyle='None',marker='^',markersize=3,linewidth=0.1)
    PLTBc,=pb.plot(B['dist'],B['c'],color=WARM,linestyle='None',marker='o',markersize=3,linewidth=0.1)
    PLTBd,=pb.plot(B['dist'],B['d'],color=WARM,linestyle='None',marker='^',markersize=3,linewidth=0.1)
    PLTCc,=pb.plot(C['dist'],C['c'],color=HOT,linestyle='None',marker='o',markersize=3,linewidth=0.1)
    PLTCd,=pb.plot(C['dist'],C['d'],color=HOT,linestyle='None',marker='^',markersize=3,linewidth=0.1)
    pb.legend([PLTAc, PLTAd,PLTBc, PLTBd,PLTCc, PLTCd], ["$\Gamma_{o/-}^{(F)}$ (1.40V)", "$\Gamma_{+/o}^{(F)}$ (1.40V)","$\Gamma_{o/-}^{(F)}$ (1.45V)", "$\Gamma_{+/o}^{(F)}$ (1.45V)","$\Gamma_{o/-}^{(F)}$ (1.50V)", "$\Gamma_{+/o}^{(F)}$ (1.50V)"],fontsize=6,borderaxespad=1.0)
    pb.text(1.45,3500,'(a)',fontsize=12)
    pb.grid(True)
    pb.ylim([0,3500])
    pb.xlim([2,4.5])
    pb.xticks(fontsize=6)
    pb.yticks(fontsize=6)
    pb.xlabel('Distance from DB (nm)',fontsize=8)
    pb.ylabel('Filling Rates (Hz)',fontsize=8,labelpad=2)
    
    #pb.tight_layout()
    #pb.savefig('./fill.eps')
    
    #EMPTYING RATES
    YA=A['a']+A['b']
    YAav=0;N=0
    for y in YA: 
      if not math.isnan(y): 
        YAav+=y; N+=1
    YAav=YAav/N
    print YAav
    YB=B['a']+B['b']
    YBav=0;N=0
    for y in YB: 
      if not math.isnan(y): 
        YBav+=y; N+=1
    YBav=YBav/N
    YC=C['a']+C['b']
    YCav=0;N=0
    for y in YC: 
      if not math.isnan(y): 
        YCav+=y; N+=1
    YCav=YCav/N
    #pb.figure(2,figsize=(3.375,3.375*3/4))
    pb.subplot(2,1,2)
    #pb.clf()
    pb.plot(np.linspace(2,4.5,300),Exp((Amp,lam),np.linspace(2,4.5,300)),'--k')
    pb.plot(np.linspace(2,4.5,300),[YAav]*300,'--',color=COOL,linewidth=1)
    pb.plot(np.linspace(2,4.5,300),[YBav]*300,'--',color=WARM,linewidth=1)
    pb.plot(np.linspace(2,4.5,300),[YCav]*300,'--',color=HOT,linewidth=1)
    PLTAa,=pb.plot(A['dist'],A['a'],color=COOL,linestyle='None',marker='s',markersize=3,linewidth=0)
    PLTAb,=pb.plot(A['dist'],A['b'],color=COOL,linestyle='None',marker='D',markersize=3,linewidth=0)
    PLTBa,=pb.plot(B['dist'],B['a'],color=WARM,linestyle='None',marker='s',markersize=3,linewidth=0)
    PLTBb,=pb.plot(B['dist'],B['b'],color=WARM,linestyle='None',marker='D',markersize=3,linewidth=0)
    PLTCa,=pb.plot(C['dist'],C['a'],color=HOT,linestyle='None',marker='s',markersize=3,linewidth=0)
    PLTCb,=pb.plot(C['dist'],C['b'],color=HOT,linestyle='None',marker='D',markersize=3,linewidth=0)
    pb.legend([PLTAa, PLTAb,PLTBa, PLTBb,PLTCa, PLTCb], ["$\Gamma_{-/o}^{(E)}$ (1.40V)", "$\Gamma_{o/+}^{(E)}$ (1.40V)","$\Gamma_{-/o}^{(E)}$ (1.45V)", "$\Gamma_{o/+}^{(E)}$ (1.45V)","$\Gamma_{-/o}^{(E)}$ (1.50V)", "$\Gamma_{o/+}^{(E)}$ (1.50V)"],fontsize=6,borderaxespad=1.0)
    pb.text(1.45,3500,'(b)',fontsize=12)
    pb.grid(True)
    pb.ylim([0,3500])
    pb.xlim([2,4.5])
    pb.xticks(fontsize=6)
    pb.yticks(fontsize=6)
    pb.xlabel('Distance from DB (nm)',fontsize=8)
    pb.ylabel('Emptying Rates (Hz)',fontsize=8,labelpad=2)
    pb.tight_layout()
    #pb.savefig('./empty.eps')
    pb.savefig('./rates.eps')
    pb.show()
  
  return
  
if __name__=='__main__':
  main()
  