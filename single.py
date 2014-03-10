#!/usr/bin/python

import sys
import ttCorr as ttC
import numpy as np
import os
import re

class get_files:
  def __init__(self):
    self.path1=sys.argv[1]                    #read path from the command line argument
    self.a=np.nan
 
  def list_files(self,a):
    for filename in os.listdir(self.path1):   #os.listdir lists files in path
#      print filename                         #print files in path  
      a.append(filename)                      
    return a

#lf=[]              # list of files
#p=get_files()      # instantiate object
#p.list_files(lf)   # get list of files to array
#print lf

def main():         # read the files and process the Gaussians
  filename=sys.argv[1]
  ttC.Analyze(filename)
       

if __name__=='__main__':
  main()
