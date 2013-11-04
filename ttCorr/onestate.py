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
    if 'A=' in var: 
      guessparams=list(guessparams)
      guessparams[0]=float(var[2:])
      guessparams=tuple(guessparams)
    elif 'mu=' in var: 
      guessparams=list(guessparams)
      guessparams[1]=float(var[3:])
      guessparams=tuple(guessparams)
    elif var=='show':
      showHist(x,n,tuple(guessparams)+(0,-1000,0,-1000))
    elif var=='y':
      good=True
    elif var!='':
      print "Use one of the accepted parameters (A,mu)"
  return guessparams

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

def oneGauss(p,x):
  A, mu = p
  sig=sigofmu(mu)
  G=(A/(sig*np.sqrt(2*np.pi)))*np.exp(-0.5*(x-mu)*(x-mu)/(sig*sig))
  return G

def oneGaussErr(p, x, n): 
  error=oneGauss(p, x) - n
  if p[0]<0:
    error=1000*error
  #minsep=1.5
  #if np.abs(p[2]-p[5])<minsep or np.abs(p[8]-p[5])<minsep or np.abs(p[8]-p[3])<minsep:
  #  error=1000*error
  return error

def oneGaussFit(x,n,guess=(3000,20),mode=None):
  params,success=leastsq(oneGaussErr,guess,args=(x,n))
  A,mu=params
  sig=sigofmu(mu)
  if mode=='vocal':
    print 'Parameters for single Gaussian fit to whole histogram:'
    print 'A: '+str('%d'%A)
    print 'sig: '+str('%.2f'%sig)
    print 'mu: '+str('%.2f'%mu)+'\n'
    
  return params
  
  