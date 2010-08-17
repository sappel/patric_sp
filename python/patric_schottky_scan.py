import os
import string
from numpy import *
from pylab import *
from scipy.fftpack import fft, fftfreq
from scipy.special.basic import erf
from obflib import smooth

# input

inputs=['15','16','32','48','47','46']
Nfiles=len(inputs)
for j in range(Nfiles):
	inputs[j]="/d/bhs01/boine/patric_schottky_bunch" + inputs[j]
inputdirs=inputs

baretune=6.0*97.3/360.0
nmode=0.0
dQb=0.0045
Seff=abs((0.5*baretune+nmode)*6.0e-3) 
smooth_factor=1
synfreq=9.1e3 # 9.1e3 # 8.2e3 
#synfreq=13.0e3   
zmin=5.5
zmax=5.5

# sc vector (x axis)

dQs_vec=zeros(Nfiles,float)
for j in range(Nfiles):
	dQs_vec[j]=dQb*j

# read idl.dat for first file only

datafile=open('%s/idl.dat' %(inputdirs[0]),'r')
lines=datafile.readlines()
datafile.close()

e_kin=float(lines[1])   
Zb=float(lines[2])      
Ab=float(lines[3])      
circum=float(lines[5])       
Nelements=float(lines[6])   
cell_length=float(lines[7])  

qe=1.6022e-19
mp=1.6605e-27
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   
freq0=beta0*clight/circum
Qs=synfreq/freq0

# arrays

Ns=int(4*4096*Nelements)
han=hanning(Ns)
s=zeros(Ns,float)
pickup_h=zeros(Ns,float)

# read first patric.dat 

datafile=open('%s/patric.dat' %(inputdirs[0]),'r')
lines=datafile.readlines()
datafile.close()

for l in xrange(Ns):
	words=string.split(lines[l])
	s[l]=float(words[0])
	pickup_h[l]=float(words[16])

Dt=(s[1]-s[0])/(beta0*clight)
freq=fftfreq(Ns,Dt)
tunevec=freq/(beta0*clight/circum)
tune_min=nmode+baretune-zmin*Qs
tune_max=nmode+baretune+zmax*Qs
jmin=tunevec.searchsorted(tune_min)
jmax=tunevec.searchsorted(tune_max)
schottky_h_sum=zeros((Nfiles,jmax-jmin),float)

# read patric.dat for each file

for j in range(Nfiles):
	datafile=open('%s/patric.dat' %(inputdirs[j]),'r')
	lines=datafile.readlines()
	datafile.close()

	for l in xrange(Ns):
		words=string.split(lines[l])
		s[l]=float(words[0])
		pickup_h[l]=float(words[16])

    # FFT
	pickup_h_fft=fft(han*pickup_h)
	schottky_h_sum[j,:]=smooth(abs(pickup_h_fft[jmin:jmax])**2,window_len=smooth_factor)
	print 'fft %s' %(j)    
	     
    # sum over +sidebands

	nmax=10  # was 10
	for u in range(nmax):
		jmin2=tunevec.searchsorted(tune_min+u+1.0)
		jmax2=jmin2+(jmax-jmin) # len(schottky_h_sum)
		schottky_h_sum[j,:]+=smooth(abs(pickup_h_fft[jmin2:jmax2])**2,window_len=smooth_factor) 
	
	schottky_h_sum[j,:]=1.0+log(schottky_h_sum[j,:])
		    
maxoffset=max(schottky_h_sum[0,:])	
Ygrid,Xgrid=meshgrid(tunevec[jmin:jmax]-baretune-nmode,dQs_vec)

Nlevels=30
cp1=contourf(Xgrid/Qs,Ygrid/Qs,schottky_h_sum,Nlevels,cmap=cm.spectral)

#for j in range(Nfiles):
#	plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,schottky_h_sum[j,:],'k')
#xlim(-zmin,zmax)
#ylim(0.0*maxoffset,1.1*maxoffset)

savefig('out.eps')
os.system('display out.eps')