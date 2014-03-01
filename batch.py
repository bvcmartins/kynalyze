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
  out=open('output.dat','w')
  lf=[]             # start the list of files
  lc=[]             # start the list of centres
  a=[]
  p=get_files()     # instantiate class get_files
  p.list_files(lf)  # object p list files in folder path given in argument line
  lf.sort()         # sort file list by file number - important for the loop
  for i in lf:
    if re.search('ANSWER',i):  # select files ANSWER
      j=sys.argv[1]+i          # concatenate path and file name
      f1=open(j,"r")           # open this file
      data=[line.split() for line in f1]  # extract lines of the file using line split
      lc=data[3],data[4],data[5]    # get the centres from the file              
      j=re.sub('_ANSWER','',i)      # get the actual data files
      j=sys.argv[1]+j          # concatenate path and file name
#      print j
      a=ttC.Analyze(j)
#      out.write(str(a[0])+' '+str(a[1])+' '+str(a[2])+' '+str(data[3])+' '+str(data[4])+' '+str(data[5])+'\n')
#def main():
#  filename=sys.argv[1]
#  ttC.Analyze2(filename)
       

if __name__=='__main__':
  main()
