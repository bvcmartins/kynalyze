#!/usr/bin/python

import sys
import numpy as np
import ttCorr as ttC

def main():
  if len(sys.argv)==1:
    print '\n/-----\\\n|ERROR|===>>> Try inputs of the form python Analyze.py foo.dat (mode)\n\-----/\n'
    sys.exit(0)
  filename=sys.argv[1]
  if len(sys.argv)==3:
    MODE=sys.argv[2] # modes are meant to produce a particular plot or particular data... ALL gets them all
  else:
    MODE=None
  if MODE in ['histogram','Irange','trace','traceCorr','corrHists',None]:
    print ''
    #print '\nRunning Analyze.py in mode "'+str(MODE)+'"\n'
  else:
    print "\nWatch it son!\nYou are only allowed modes in:\n['histogram','Irange','trace',traceCorr,corrHists,None]\n"
    sys.exit(0)
  ttC.Analyze(filename)
  return

if __name__=='__main__':
  main()
