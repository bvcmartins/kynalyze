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
from operator import itemgetter
import re

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

def smooth(TR):
  window_len=11
  window='hanning'
  s=np.r_[TR.n[window_len-1:0:-1],TR.n,TR.n[-1:-window_len:-1]]
  w=eval('np.'+window+'(window_len)')
  TR.y=np.convolve(w/w.sum(),s,mode='valid')
  TR.y=TR.y[5:-5]
#  pb.plot(TR.x,TR.y,color=col1,linewidth=1)
#  pb.show() 

def peak_finder(TR):
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]

  for i in np.where(pk==1)[0]:
    TR.peak.append(TR.x[int(i)])
    TR.amp.append(TR.y[int(i)])

#  TR.amp=zip(TR.amp,TR.peak)
#  TR.amp=sorted(TR.amp, reverse=True)

#  npeak=0
#  n=0
#  tmp=np.zeros((10,2))
#  for i in pk:
#    if i==1:
##      TR.peak.append(TR.x[n])
##      TR.amp.append(TR.y[n])
#      tmp[npeak][0]=TR.x[n]    # peak position
#      tmp[npeak][1]=TR.y[n]    # amplitude
#      npeak+=1   
#    n+=1
#  tmp=sorted(tmp,key=itemgetter(1),reverse=True)
#  tmp=tmp[0:3]
#  tmp=sorted(tmp,key=itemgetter(0))
#  TR.peak=tmp[0][0],tmp[1][0],tmp[2][0]
#  TR.amp=tmp[0][1],tmp[1][1],tmp[2][1]

#  print TR.peak
#  print TR.amp
#  pb.plot(TR.x,TR.y,color=col1,linewidth=1)
  pb.bar(TR.x,TR.n,float(TR.Imax)/TR.numbin,color='0.8',linewidth=0.4,align='center')
#  pb.show()
#  raise sys.exit()

def get_guessparams(TR):

#  PoNEG=0.33
#  PoNEU=0.33
#  PoPOS=0.33
#  mu2=TR.x[list(TR.n).index(max(TR.n))]
#  mu1=0.9*mu2
#  mu3=1.1*mu2  
#  guessparams=(PoNEG*2*TR.Fs*TR.Imax/TR.numbin, mu1, PoNEU*2*TR.Fs*TR.Imax/TR.numbin, mu2, PoPOS*2*TR.Fs*TR.Imax/TR.numbin, mu3)
  guessparams=(TR.amp[0],TR.peak[0],TR.amp[1],TR.peak[1],TR.amp[2],TR.peak[2])
  return guessparams

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

def threeGaussFit(TR):
  print len(TR.peak)

  if len(TR.peak)>0:
    if len(TR.peak)==1:
      sig=[ttC.sigofmu(TR.peak[0])]
      popt,pcov=curve_fit(oneGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0]])
      TR.popt.append(popt[0])
      TR.popt.append(popt[1])
      TR.popt.append(popt[2])
#      pb.plot(TR.x,oneGauss2(TR.x,popt[0],popt[1],popt[2]),color=col4,linewidth=1)
#      pb.show()
      TR.n=TR.n-oneGauss2(TR.x,popt[0],popt[1],popt[2])

    elif len(TR.peak)==2:
      sig=[ttC.sigofmu(TR.peak[0]),ttC.sigofmu(TR.peak[1])]
      popt,pcov=curve_fit(twoGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1]])

      for i in range(6):
        TR.popt.append(popt[i])
  
#      pb.plot(TR.x,twoGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),color=col4,linewidth=1)
#      pb.show()
      TR.n=TR.n-twoGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])

    else:
      sig=[ttC.sigofmu(TR.peak[0]),ttC.sigofmu(TR.peak[1]),ttC.sigofmu(TR.peak[2])]
      popt,pcov=curve_fit(threeGauss2,TR.x,TR.n,p0=[TR.peak[0],TR.amp[0],sig[0],TR.peak[1],TR.amp[1],sig[1],TR.peak[2],TR.amp[2],sig[2]])
      for i in range(9):
        TR.popt.append(popt[i])
  
#      pb.plot(TR.x,threeGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8]),color=col4,linewidth=1)
#      pb.show()
      TR.n=TR.n-threeGauss2(TR.x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8])

def hidden_peak(TR,filename):
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
      sig=[ttC.sigofmu(peak[0])]
      popt,pcov=curve_fit(oneGauss2,x,n,p0=[peak[0],amp[0],sig[0]])
      TR.popt.append(popt[0])
      TR.popt.append(popt[1])
      TR.popt.append(popt[2])

    elif len(peak)==2:
      sig=[ttC.sigofmu(peak[0]),ttC.sigofmu(peak[1])]
      popt,pcov=curve_fit(twoGauss2,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1]])

      for i in range(6):
        TR.popt.append(popt[i])
    
    else:
      sig=[ttC.sigofmu(peak[0]),ttC.sigofmu(peak[1]),ttC.sigofmu(peak[2])]
      popt,pcov=curve_fit(threeGauss2,x,n,p0=[peak[0],amp[0],sig[0],peak[1],amp[1],sig[1],peak[2],amp[2],sig[2]])
 
      for i in range(9):
        TR.popt.append(popt[i])
#
###

  if len(TR.popt)==3:
     pb.plot(x,oneGauss2(x,TR.popt[0],TR.popt[1],TR.popt[2]),color=col3,linewidth=1)
  elif len(TR.popt)==6:
     pb.plot(x,twoGauss2(x,TR.popt[0],TR.popt[1],TR.popt[2],TR.popt[3],TR.popt[4],TR.popt[5]),color=col3,linewidth=1)
  else:
     pb.plot(x,threeGauss2(x,TR.popt[0],TR.popt[1],TR.popt[2],TR.popt[3],TR.popt[4],TR.popt[5],TR.popt[6],TR.popt[7],TR.popt[8]),color=col3,linewidth=1)

###
  i=re.search('[0-9]+',filename)
  pb.savefig('figures/'+i.group()+'.png')
  pb.clf()

  return TR.popt


#  TR.popt+=[np.float64(0.0),np.float64(0.0),np.float64(1.0)]  # this is provisory. I want to find a better solution.
#  TR.popt+=[np.float64(0.0),np.float64(0.0),np.float64(1.0)]
#
##    for i in range(3):
##      if popt[3*i] > 0 and popt[3*i+1]>TR.sens:
##        np.append(TR.popt,popt[3*i])
##        np.append(TR.popt,popt[3*i+1])
##        np.append(TR.popt,popt[3*i+2])
##
#  print TR.popt
#
#  pb.plot(x,threeGauss2(x,TR.popt[0],TR.popt[1],TR.popt[2],TR.popt[3],TR.popt[4],TR.popt[5],TR.popt[6],TR.popt[7],TR.popt[8]),color=col3,linewidth=1)
#  i=re.search('[0-9]+',filename)
#  pb.savefig('figures/'+i.group()+'.png')
#  pb.clf()
#
#  return TR.popt
#
