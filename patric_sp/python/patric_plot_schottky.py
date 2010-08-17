import os
import string
from numpy import *
from pylab import *
from scipy.fftpack import fft, fftfreq
from scipy.special.basic import erf
from pickle import Pickler, Unpickler
#from obflib import smooth

# input

inputdir="/u/boine/d/patric_decoh_bunch1"
#inputdir="/d/bhs01/boine/patric_schottky_dc5"
#inputdir="/u/boine/d/patric_schottky_b_new3"

baretune=6.0*97.3/360.0
nmode=0.0
Seff=abs((0.5*baretune+nmode)*6.0e-3) 
smooth_factor=1
Usc=5.0
Qbsc=0.0056*1.0  
#Qbsc=0.0085*2.0           
Qc=0.0045*0.0  # 0.008*0.379
synfreq=8.2e3 # 9.1e3 # 9.1e3 # 8.2e3 
#synfreq=13.0e3  
Npmax=0 #4*4096 
zmin=6.0
zmax=6.0

Usc=Qbsc/Seff
print "Usc:",Usc


#---

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

qe=1.6022e-19
mp=1.6605e-27
clight=2.99e8
gamma0 = 1.0 + (e_kin*1e6*qe)/(mp*clight*clight) 
beta0  = sqrt((gamma0*gamma0-1.0)/(gamma0*gamma0))   
freq0=beta0*clight/circum
Qs=synfreq/freq0
print "Qs:", Qs
print "dQsc/Qs:", Qbsc/Qs

l_syn=beta0*clight/synfreq

# read patric.dat

datafile=open('%s/patric.dat' %(inputdir),'r')
lines=datafile.readlines()
datafile.close()

Np=len(lines)
print 'cells:',Np/Nelements

if Npmax > 0 and Npmax < Np/Nelements :
  Np=int(Npmax*Nelements)

s=zeros(Np,float)
xz=zeros(Np,float)
xz2=zeros(Np,float)
xrms=zeros(Np,float)
pickup_h=zeros(Np,float)
offset_h=zeros(Np,float)
zx=zeros(Np,float)

for j in xrange(Np):
     words=string.split(lines[j])
     s[j]=float(words[0])
     xrms[j]=float(words[5])
     xz2[j]=float(words[8])
     xz[j]=float(words[9])
     offset_h[j]=float(words[13])
     pickup_h[j]=float(words[16])

# FFT

han=hanning(Np)
offset_h_fft=fft(han*offset_h)
pickup_h_fft=fft(han*pickup_h)
xrms_fft=fft(han*xrms)
xz_fft=fft(han*xz)
Dt=(s[1]-s[0])/(beta0*clight)
freq=fftfreq(Np,Dt)
tunevec=freq/(beta0*clight/circum)

"""
fpick=open('pick.dat','w')
p1=Pickler(fpick)
p1.dump(pickup_h_fft)
fpick.close()

fpick=open('pick.dat','r')
u1=Unpickler(fpick)
pickup_h_fft2=u1.load()
fpick.close()
"""

# analytic Schottky

def c_erf(z):
	return exp(-z**2)*(1.0-erf(-z*1.0j))

def schottky_dc(tune,baretune,Seff,Usc):
	z=(tune-baretune)/Seff+Usc
	return exp(-0.5*z**2)/abs(1.0+sqrt(0.5*pi)*Usc*1.0j*c_erf(z/sqrt(2.0)))**2

def schottky_burov(tune,baretune,Seff, Usc):
	Usc2=Usc
	lambda1=sqrt(0.5*pi)*Usc2**2*Seff*exp(-0.5*Usc2**2-1.0)
	dQ2=Seff/Usc2
	return lambda1/((tune-baretune-dQ2)**2+lambda1**2)

def barrier_ab(dQsc, Qs, m):
	dQm=0.0
	if m > 0.0 : 
		dQm=-(dQsc+Qc)/2.0+sqrt((dQsc/2.0-Qc/2.0)**2+(m*Qs)**2)
	else:
		dQm=-(dQsc+Qc)/2.0-sqrt((dQsc/2.0-Qc/2.0)**2+(m*Qs)**2)
	return dQm/Qs

def long_WB(tune, Seff):
	n=size(tune)
	fd=zeros(n,float)
	for i in range(n):
		fd[i]=0.0 
		if 1.0-tune[i]**2/Seff**2 > 0.0:
			fd[i]=(1.0-tune[i]**2/Seff**2)
	return fd
	

# plot

tune_min=nmode+baretune-zmin*Qs
tune_max=nmode+baretune+zmax*Qs
jmin=tunevec.searchsorted(tune_min)
jmax=tunevec.searchsorted(tune_max)
schottky_h_sum=smooth(abs(pickup_h_fft[jmin:jmax])**2,window_len=smooth_factor)
#schottky_h_sum2=smooth(abs(pickup_h_fft2[jmin:jmax])**2,window_len=smooth_factor)

# sum over +sidebands

nmax=6  # was 10
for j in range(nmax):
	  jmin2=tunevec.searchsorted(tune_min+j+1.0)
	  #jmax2=tunevec.searchsorted(tune_max+j+1.0)
	  jmax2=jmin2+len(schottky_h_sum)
	  schottky_h_sum+=smooth(abs(pickup_h_fft[jmin2:jmax2])**2,window_len=smooth_factor)
	  #schottky_h_sum2+=smooth(abs(pickup_h_fft2[jmin2:jmax2])**2,window_len=smooth_factor)
	
	 
#maxoffset=max(smooth(abs(pickup_h_fft[jmin:jmax])**2,window_len=smooth_factor))
maxoffset=max(schottky_h_sum)
#maxoffset=max(smooth(abs(offset_h_fft[jmin:jmax])**2,window_len=smooth_factor))
Seff2=abs((0.5*baretune+1.0)*6.0e-3) 
schott_ana=schottky_dc(tunevec[jmin:jmax],baretune,Seff,Usc)
norm0=0.0
for l in range(jmax-jmin):
	norm0+=schott_ana[l]
schott_ana2=schottky_burov(tunevec[jmin:jmax],baretune,Seff,Usc)
fWB=long_WB(tunevec[jmin:jmax]-baretune-nmode+Qbsc,sqrt(5.0)*Seff2)
maxana2=max(schott_ana2)

for j in range(nmax):
	Seff2=abs((0.5*baretune+j+1.0)*6.0e-3) 
	norm2=0.0
	schott_ana_tmp=schottky_dc(tunevec[jmin:jmax],baretune,Seff2,Usc*Seff/Seff2)
	for l in range(jmax-jmin):
		norm2+=schott_ana_tmp[l]
	schott_ana+=norm0/norm2*schott_ana_tmp

maxana=max(schott_ana**2)

print jmin, jmax

"""
for j in xrange(schottky_h_sum.shape[0]):
     if (schottky_h_sum[j] > 0.02):
	    print "tune_m, P_m:",(tunevec[jmin+j]-baretune-nmode)/Qs, schottky_h_sum[j] 
"""
#semilogy

subplot(211)

plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,schottky_h_sum,'k',label=r'$\Delta Q_{sc}/Q_s=2$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,abs(offset_h_fft[jmin:jmax]**2),'k',label='xz')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Seff,smooth(abs(offset_h_fft[jmin:jmax])**2,window_len=smooth_factor),'k',label=r'$\rm{simulation}$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,schottky_h_sum2,'r',label=r'$\Delta Q_{sc}/Q_s=0$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,smooth(abs(pickup_h_fft_sum)**2,window_len=smooth_factor),'k',label='WB')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,smooth(abs(pickup_h_fft2[jmin:jmax])**2,window_len=smooth_factor),'r',label='KV')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,0.8*maxoffset*fWB,'k:',label='$\mathrm{dc\;beam}$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs-0.0*Qc/Qs,1.0*maxoffset/maxana*schott_ana**2,'g:',label=r'$S(f)$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Seff,maxoffset*exp(-0.5*(tunevec[jmin:jmax]-baretune-nmode)**2/(Seff)**2),'r',label=r'$S_0(f)$')
#plot((tunevec[jmin:jmax]-baretune-nmode)/Qs,maxoffset*exp(-0.5*(tunevec[jmin:jmax]-baretune-nmode)**2/(Seff)**2),'r',label=r'$\rm{Gauss:} U_{sc}=0$')
#plot(tunevec[jmin:jmax],maxoffset/maxana2*schott_ana2,'g')

dQsc=-Seff*Usc
axvline(-Qbsc/Qs,linewidth=2,color='r',linestyle=':',label=r'$Q_0-\Delta Q_{sc}$')
#axvline(dQsc/Seff,linewidth=2,color='k',linestyle=':',label=r'$Q_0-\Delta Q_{sc}$')
axvline(-Qc/Qs,color='k',linestyle=':',label=r'$Q_0$')
#axvline(-Qc/Seff,color='k',linestyle='--',label=r'$Q_0$')


axvline(barrier_ab(Qbsc,Qs,1),color='b',linestyle='--',label=r'$\Delta Q_{k=\pm 1}$')
axvline(barrier_ab(Qbsc,Qs,-1),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,2),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-2),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,3),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-3),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,4),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-4),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,5),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-5),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,6),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-6),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,7),color='b',linestyle='--')
axvline(barrier_ab(Qbsc,Qs,-7),color='b',linestyle='--')

"""
axvline(barrier_ab(0.0,Qs,1),color='k',linestyle='--',label=r'$\Delta Q_{k=\pm 1}$')
axvline(barrier_ab(0.0,Qs,-1),color='k',linestyle='--')
"""

xlim(-zmin,zmax)
ylim(0.0*maxoffset,1.1*maxoffset)
xlabel(r'$(Q_f-Q_0)/Q_s$')
#xlabel(r'$(Q_f-Q_0)/\delta Q_{\xi}$')
ylabel('S [arb. units]')

#legend(loc=2)  
subplots_adjust(left=0.18,bottom=0.12)
#subplots_adjust(top=0.95,hspace=0.35,bottom=0.075,left=0.15)



subplot(212)

s_red=s[0:Np-1:50]
xz2_red=xz2[0:Np-1:50]
xz_red=xz[0:Np-1:50]
#pickup_h_red=pickup_h[0:Np-1:40]
#offset_h_red=offset_h[0:Np-1:40]
#plot(s_red/cell_length,pickup_h_red,'k')
#plot(s_red/cell_length,offset_h_red,'r')
#plot(s_red/(cell_length*12.0),xz2_red,'r')
#plot(s_red/(cell_length*12.0),xz_red,'b')
plot(s_red/l_syn,xz2_red,'r')
plot(s_red/l_syn,xz_red,'b')
xlabel('turns')
ylabel('offset')
subplots_adjust(left=0.18,bottom=0.12)


show()
#savefig('out.eps')
#os.system('display out.eps')
    