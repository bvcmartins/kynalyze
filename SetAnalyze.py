import sys
import os
import numpy as np
import pylab as pb
import scipy.fftpack
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import ttCorr as ttC
import Analyze
import pickle

col1='#FF7A00'
col2='#03899C'
col3='#D2006B'
col4='#619A00'

numbin=300
Imax=100
Fs=1e4

def main():
  if len(sys.argv)==1:
    print '\n/-----\\\n|ERROR|===>>> Try inputs of the form python Analyze.py foo.p where foo.p is a pickle... not literally.\n\-----/\n'
    sys.exit(0)
  picklename=sys.argv[1]
  SET=[]
  x=np.array([0.0000000020128,0.000000002108307,0.000000002205033,0.000000002304373,0.000000002404483,0.000000002506767,0.000000002609444,0.000000002714023,0.000000002819564,0.000000002925106,0.000000003032267,0.000000003139197,0.000000003247665,0.000000003356696,0.000000003465344,0.00000000357534,0.000000003684857,0.000000003795656,0.000000003905897,0.000000004017367,0.000000004129093,0.000000004240215,0.000000004352506,0.000000004464105])
  x=[i*10**9 for i in x]
  basename='./20130506_06/A'
  filenumbers=range(101)[40:71]#[6:20]
  filenumbers=[('%03d'%filenumber) for filenumber in filenumbers]
  filenames=[basename+filenumber+'.dat' for filenumber in filenumbers]
  #GET RID OF THIS NEXT!
  #filenames=filenames[12:15]
  for filename in filenames: print filename
  print ''
  for filename in filenames:
    
    success=False
    while success==False:
      TR=Analyze.Analyze(filename)
      succ=False
      while succ==False:
        IN=raw_input("y to accept, or 'skip' or 'repeat': ")
        if IN=='y':
          succ=True
          SET.append(TR)
          success=True
        elif IN=='skip':
          succ=True
          success=True
        elif IN=='repeat':
          succ=True
        elif IN!='':
          print "What!?"
        
  pickle.dump(SET,open( picklename, 'w'))
  
  pb.show()
  
  
  
  #PTa,=pb.plot(x,alist,col1,marker='o',markersize=12,markerfacecolor=col1,linewidth=2,linestyle='-')
  #PTb,=pb.plot(x,blist,col2,marker='^',markersize=12,markerfacecolor=col2,linewidth=2,linestyle='--')
  #PTc,=pb.plot(x,clist,col3,marker='v',markersize=12,markerfacecolor=col3,linewidth=2,linestyle=':')
  #PTd,=pb.plot(x,dlist,col4,marker='<',markersize=12,markerfacecolor=col4,linewidth=2,linestyle='-.')
  #pb.grid(True)
  #pb.xticks(fontsize=14)
  #pb.yticks(fontsize=14)
  #pb.xlabel('distance from DB (nm)',fontsize=18)
  #pb.ylabel('rates (Hz)',fontsize=18)
  #pb.legend([PTa,PTb,PTc,PTd],["a", "b","c","d"],fontsize=18)
  #pb.show()

if __name__ == '__main__':
  main()



