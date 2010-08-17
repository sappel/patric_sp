import os
import string
from numpy import *
from Scientific.IO.ArrayIO import *
from pylab import *
from TransLandau import w_gauss, parabolic, elliptic

outputdir='/d/bhs01/boine/track_sc_dc_el_n10Q2_fine_new/'
# x parameter:
Xpara='Rs'
Ypara='Zimage'

# count data sets in dir:
files=os.listdir(outputdir)
Njobs=len(files)
print "N jobs in dir:", Njobs

Xvec=zeros(Njobs,float32)
Yvec=zeros(Njobs,float32)
Zvec=zeros(Njobs,float32)

#store parameters and data in vecs:
for jobid in xrange(1,Njobs+1):
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
  #datavec=readFloatArray('%srun%s/patric.dat' %(outputdir,jobid))
  #pickup_h=datavec[:,13]
  datafile=open('%srun%s/patric.dat' %(outputdir,jobid),'r')
  lines=datafile.readlines()
  datafile.close()
  pickup_h=zeros(len(lines),float)
  for j in xrange(len(lines)):
      words=string.split(lines[j])
      pickup_h[j]=float(words[13])   # 1: emittance_h, 13: offset_h
  Zvec[jobid-1]=max(pickup_h)-pickup_h[0]
  #if Zvec[jobid-1] >= 0.002 and Zvec[jobid-1] < 0.02 :
  #	  Zvec[jobid-1] = 0.02 
  print "jobid:", jobid
  
# coordinate vectors:

Xcoord=[]
Ycoord=[]

# sort copies of coordinate vecs

Xsort=Xvec.copy()
Xsort.sort()
Ysort=Yvec.copy()
Ysort.sort()

tmp=Xsort[0]
Xcoord.append(tmp)
for x in Xsort:
   if x != tmp:
      Xcoord.append(x)	   
   tmp=x
tmp=Ysort[0]
Ycoord.append(tmp)
for x in Ysort:
   if x != tmp:
      Ycoord.append(x)	   
   tmp=x

# coordinate grid   
Nx=len(Xcoord)
Ny=len(Ycoord)

# function array to be plotted
Zfunc=zeros((Nx,Ny),float32)
for j in xrange(0,Njobs):
    l=Xcoord.index(Xvec[j])
    k=Ycoord.index(Yvec[j])
    Zfunc[l,k]=Zvec[j] 

print Zfunc   

# convert to numpy arrays for plotting:
Xcoord=array(Xcoord)
Ycoord=array(Ycoord)
Ygrid,Xgrid=meshgrid(Ycoord,Xcoord)

# read analytic stability boundary from files

ratef=0.001j  # damping rate

zfreq=arange(-20.0,20.0,0.1,float)   # argument z runs from -10 to 20 with step 0.25
zfreq2=arange(-1.1,1.1,0.0001,float)
Npoints=size(zfreq) 
Npoints2=size(zfreq2) 
uv1=zeros(Npoints,complex)
uv2=zeros(Npoints,complex)
uv3=zeros(Npoints,complex)
uv4=zeros(Npoints2,complex)
uv5=zeros(Npoints,complex)

for j in range(Npoints):	
   uv1[j]=1.0/w_gauss(zfreq[j]+ratef)
   uv5[j]=1.0/w_gauss(zfreq[j]+0.1j)
   uv2[j]=sqrt(5.0)*parabolic(zfreq[j]) 
   #uv3[j]=-1.0j/gauss_octupole(zfreq[j]-ratef)
   uv5[j]=complex(j*2.0/Npoints,-25.0*j*2.0/Npoints)
for j in range(len(zfreq2)):
   uv4[j]=sqrt(5.0)*elliptic(zfreq2[j])
  
   
#------plot-----------


print "start plot"

params = {'backend': 'ps',
	  'font.size:': 20,
          'axes.labelsize': 24,
	  'axes.linewidth':2,
         'text.fontsize': 20,
	 'legend.fontsize': 20, 
         'xtick.labelsize': 24,
	 'xtick.major.size': 6, 
	 'xtick.minor.size': 4,
	 'xtick.major.pad': 8, 
	 'xtick.minor.pad': 8,
         'ytick.labelsize': 24,
	 'ytick.major.size': 6, 
	 'ytick.minor.size': 4,
	 'ytick.major.pad': 8, 
	 'ytick.minor.pad': 8,
	 'lines.linewidth': 3,
	 'figure.figsize': (10,7),
	 'text.usetex': False
	 }

	 
rcParams.update(params)

subplot(111)
Nlevels=40
factorx=1.4e7*1.0e-6  #1.4e7*1.0e-6  
factory=1.4e7*1.0e-6  #5.6e6*1.0e-6 
cp1=contourf(1.0e-6*Xgrid,1.0e-6*Ygrid,Zfunc,Nlevels,cmap=cm.bone) 
#p1=plot(factorx*uv1.real,factorx*uv1.imag-4.0*factorx,'r',label='Stability boundary: Gauss')
p2=plot(factorx*uv4.real,factorx*uv4.imag-4.0*factorx,'r',label='dc stability boundary' )
#p3=plot(uv5.real,uv5.imag,'r--',label='image current' )
#axvline(0.2,linewidth=2,color='r',linestyle=':')
#plot(Xgrid/8.0e6,Zfunc,'bs')
plot([0.3],[-10.0],'ro',markersize=10,label='SIS 18: dc, 200 MeV/u')
cb1 = colorbar(cp1, shrink=0.8, extend='both')
#text(23.0,-85.0,'final dipole amplitude [arb. units]',rotation='vertical', size=20)
text(2.3,-45.0,'final dipole amplitude [arb. units]',rotation='vertical', size=20)
xlabel(r'$Re(Z_x)\/\/\rm{M}\Omega/m$')
ylabel(r'$Im(Z_x)\/\/\rm{M}\Omega/m$')
ylim(-50.0,0.0)  # -90.0
xlim(0.0,1.75)  # 17.5 
#legend(loc=1)
subplots_adjust(bottom=0.12,right=0.75)

#axis([0.0e7, 1.0e7, -5.0e7, 0.0e7])
savefig('out.eps')
#os.system('display out.eps &')


