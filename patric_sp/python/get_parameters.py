import os
import string
from numpy import *

outputdir='/d/bhs01/boine/track_sc_dc_el_n10Q2_fine_new/'
print outputdir
# x parameter:
Xpara='Rs'
Ypara='Zimage'

# count data sets in dir:
files=os.listdir(outputdir)
Njobs=len(files)
print "N jobs in dir:", Njobs

Xvec=zeros(Njobs,float32)
Yvec=zeros(Njobs,float32)

#store parameters and data in vecs:
for jobid in range(1,Njobs+1):
  conffile=open('%srun%s/patric.cfg' %(outputdir,jobid),'r')
  lines=conffile.readlines()
  conffile.close()   
  for line in lines:
     if string.find(line,Xpara) != -1:	     
       words=string.split(line)
       Xvec[jobid-1]=float(words[1])
     if string.find(line,Ypara) != -1:
       words=string.split(line)
       Yvec[jobid-1]=float(words[1]) 
  print "jobid:",jobid,"  ",Xpara,":",Xvec[jobid-1],"   ",Ypara,":",Yvec[jobid-1]
