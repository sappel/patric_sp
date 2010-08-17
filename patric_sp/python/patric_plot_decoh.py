import os
import string
import struct
from numpy import *
from pylab import *
from scipy.fftpack import fft, fftfreq
from scipy.special.basic import erf
from pickle import Pickler, Unpickler
#from obflib import smooth

#---
import numpy 

def smooth(x,window_len=10,window='hanning'):
	
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
 
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
  
	if window_len<3:
		return x
 
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
  
	s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
     #print(len(s))
	if window == 'flat': #moving average
		w=ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')
 
	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]

#----


# input

inputdir="/u/boine/d/patric_decoh_bunch1"

baretune=6.0*97.3/360.0
nmode=0.0
Seff=abs((0.5*baretune+nmode)*6.0e-3) 
smooth_factor=1
#Qbsc=0.0045*5.0  
Qbsc=0.008*5.0           
Qc=0.008*0.0  # 0.008*0.379
synfreq=8.2e3 
#synfreq=13.0e3  

#---- 

# read idl.dat

datafile=open('%s/idl.dat' %(inputdir),'r')
lines=datafile.readlines()
datafile.close()

numprocs=float(lines[0])
e_kin=float(lines[1])   
Zb=float(lines[2])      
Ab=float(lines[3])      
current=float(lines[4]) 
circum=float(lines[5])       
Nelements=float(lines[6])   
cell_length=float(lines[7])  
NZ=int(lines[12]) 

qe=1.6022e-19
mp=1.6605e-27
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   
freq0=beta0*clight/circum
Qs=synfreq/freq0
print "Qs:", Qs
print "dQsc/Qs:", Qbsc/Qs


def barrier_ab(dQsc, Qs, m):
	dQm=0.0	
	if m > 0.0 : 
		dQm=-dQsc/2.0+sqrt((dQsc/2.0)**2+(m*Qs)**2) 
	else:
		dQm=-dQsc/2.0-sqrt((dQsc/2.0)**2+(m*Qs)**2)
	return dQm/Qs

Qk=baretune+0.96*barrier_ab(Qbsc,Qs,3)*Qs


# read vector dipole_x.dat

Nt=2
dipole_x=fromfile(file='%s/dipole_x.dat' %(inputdir),dtype='f')
rho_z=fromfile(file='%s/rho_z.dat' %(inputdir),dtype='f')

Nmax=size(dipole_x)/NZ
print "Nmax:",Nmax

s_vec=zeros(Nmax)
for j in xrange(Nmax):
	s_vec[j]=120.0*cell_length*j

z_vec=zeros(NZ)
for j in xrange(NZ):
	z_vec[j]=circum/NZ*j

# norm

norm_vec=zeros(Nmax)
norm0=0.0
for j in xrange(NZ):
	norm0+=(dipole_x[j])**2

for l in xrange(Nmax):
     for j in xrange(NZ):
		norm_vec[l]+=dipole_x[j]*dipole_x[l*NZ+j]/norm0

for l in xrange(Nmax):
	if norm_vec[l] > 0.25: print "lmax:", l,"norm:",norm_vec[l]

offset_x=zeros(NZ*Nmax)
for j in range(Nmax*NZ):
	if rho_z[j] > 0.0:
		offset_x[j]=dipole_x[j]/rho_z[j]

#plot(z_vec,dipole_x[0:NZ],'k')
#plot(z_vec,1.5*dipole_x[Nt*NZ:(Nt+1)*NZ],'r')
plot(z_vec,offset_x[0:NZ],'k')
plot(z_vec,offset_x[Nt*NZ:(Nt+1)*NZ],'r')
xlim(0.0,circum)

#plot(s_vec/circum,norm_vec,'k')
#ylim(-2.0,2.0)

#savefig('out.eps')
#os.system('display out.eps')

show()

