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
#  pb.plot(TR.x,TR.n,color=col3,linewidth=1)
#  pb.show() 

def peak_finder(TR):
  pk=np.r_[1, TR.y[1:] > TR.y[:-1]] & np.r_[TR.y[:-1] > TR.y[1:],1]
  npeak=0
  n=0
  tmp=np.zeros((10,2))
  for i in pk:
    if i==1:
#      TR.peak.append(TR.x[n])
#      TR.amp.append(TR.y[n])
      tmp[npeak][0]=TR.x[n]    # peak position
      tmp[npeak][1]=TR.y[n]    # amplitude
      npeak+=1   
    n+=1
  tmp=sorted(tmp,key=itemgetter(1),reverse=True)
  tmp=tmp[0:3]
  tmp=sorted(tmp,key=itemgetter(0))
  TR.peak=tmp[0][0],tmp[1][0],tmp[2][0]
  TR.amp=tmp[0][1],tmp[1][1],tmp[2][1]
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

def threeGauss2(x,A1,mu1,sig1,A2,mu2,sig2,A3,mu3,sig3):
#  A1, mu1, sig1, A2, mu2, sig2, A3, mu3, sig3= p
#  sig1=ttC.sigofmu(mu1)
#  sig2=ttC.sigofmu(mu2)
#  sig3=ttC.sigofmu(mu3)
  G1=(A1/(sig1*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu1)*(x-mu1)/(sig1*sig1))
  G2=(A2/(sig2*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu2)*(x-mu2)/(sig2*sig2))
  G3=(A3/(sig3*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu3)*(x-mu3)/(sig3*sig3))
  return G1+G2+G3

def threeGaussFit(TR,filename):
#  guess=get_guessparams(TR)
  x=TR.x
  n=TR.n
#  params,success=leastsq(threeGaussErr,guess,args=(x,n))
#  pb.plot(TR.x,threeGauss(params,TR.x),color=col2,linewidth=1)

  mu1=TR.peak[0]
  mu2=TR.peak[1]
  mu3=TR.peak[2]
  sig1=ttC.sigofmu(mu1)
  sig2=ttC.sigofmu(mu2)
  sig3=ttC.sigofmu(mu3)
  popt,pcov=curve_fit(threeGauss2,x,n,p0=[TR.amp[0],TR.peak[0],sig1,TR.amp[1],TR.peak[1],sig2,TR.amp[2],TR.peak[2],sig3])
#  print popt
  pb.plot(x,threeGauss2(x,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6],popt[7],popt[8]),color=col3,linewidth=1)
#  pb.show()
  i=re.search('[0-9]+',filename)
  pb.savefig('figures/'+i.group()+'.png')
  pb.clf()
  return mu1,mu2,mu3
