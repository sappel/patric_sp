"""
A python script to launch PATRIC jobs 
"""

from math import *
from odict import * 
from socket import gethostname
import os

outputdir='/d/bhs01/boine/patric_schottky_dc1/'   # should exist  
runid=1                   # run identification number (only out.dat)

#-----MPI options------------------

num_procs=2
machines=['lxir033']

#----------------------------------

patric_dict=OrderedDict()

patric_dict['NPIC']=100000 
patric_dict['NX']=128
patric_dict['NY']=128
patric_dict['NZ']=256
patric_dict['NZ_bunch']=64
patric_dict['cells']=4*4096 
patric_dict['e_kin']=200.0
patric_dict['Z']=28.0
patric_dict['A']=238.0
patric_dict['current']=0.55*1.0    #        0.438*2.0 
patric_dict['piperadius']=0.05
patric_dict['circum']=108.0 
patric_dict['gamma_t']=5.4
patric_dict['CF_advance_h']=97.3*pi/180.0
patric_dict['CF_advance_v']=97.3*pi/180.0 # 129.3*pi/180.0
patric_dict['CF_R']=0.0
patric_dict['CF_length']=18.0
patric_dict['NCF']=16
patric_dict['koct']=8.0
patric_dict['dQxm']=-0.0054*1.0
patric_dict['dQym']=-0.0054*1.0
patric_dict['dqx_detune']=-0.0025/6.0   # rms
patric_dict['dqy_detune']=-0.0025/6.0   # rms
patric_dict['pic_subset']=10000
patric_dict['init_pic_xy']=1  # 0 (WB), 1 (KV), 2 (SG), 3 (GS)   # was 0 !!!!
patric_dict['init_pic_z']=2   # 1 (elliptic bunch), 0 (elliptic, coast), 2 (Gauss, coast)   # was 4 !!!!
patric_dict['momentum_spread']=6.0e-3   # 6.0e-3 rms
patric_dict['rms_emittance_x0']=5.0
patric_dict['rms_emittance_y0']=5.0
patric_dict['mismatch_x']=1.0
patric_dict['mismatch_y']=1.0
patric_dict['offcenter']=0.0  # !!! meters
patric_dict['bunchfactor']=1.0  # !!!
patric_dict['dqci']=0.0
patric_dict['dqcr']=-0.0
patric_dict['Rs']=0.0 # 8.5e6
patric_dict['nres']=10.0
patric_dict['Qs']=2.0 
patric_dict['leit']=1.0e5 
patric_dict['Zimage']=-7.0e7 # -3.0e7/2.0 
patric_dict['madx_input_file']=0
patric_dict['space_charge']=0   #   0 (off), 1 (self-consistent), 2 (linear), 3 (nonlinear)
patric_dict['imp_kick']=0       #   0 (off), 1 (on) 
patric_dict['sliced']=0
patric_dict['cavity']=0   # 0 (off), 1 (rf), 2 (barrier)
patric_dict['octupole_kick']=0
patric_dict['ampdetun_kick']=0
patric_dict['chroma']=1  
patric_dict['bc_end']=1    #  0 (open) , 1 (periodic)
patric_dict['print_cell']=120  # 120
patric_dict['footprint']=0
patric_dict['btf']=0   # was 0 !!!!!!
patric_dict['btf_harmonic']=0
patric_dict['ausgabe']=outputdir  

#-------------------run-------------------------------------

# create new PATRIC configuration file	
outfile=open('%s/patric.cfg' % (outputdir),'w')
# write PATRIC configuration file in output dir:
for x,y in patric_dict.iteritems():
   outfile.write('%s:  %s\n' %(x,y))
outfile.close()
# mpi machine configuation file:
machinefile=open('machines.tmp','w')
machinefile.write('%s:4\n' %(machines[0]))
#machinefile.write('%s:4\n' %(machines[1]))
machinefile.close()
# start PATRIC/mpi job:
thishost=gethostname()
if thishost == machines[0]:  
   cmd="nohup mpirun2 -np %d -machinefile machines.tmp ../bin/patric %s/patric > %s/out.%s.dat &" %(num_procs,outputdir,outputdir,runid)
else:
   cmd="nohup mpirun2 -nolocal -np %d -machinefile machines.tmp ../bin/patric %s/patric > %s/out.%s.dat &" %(num_procs,outputdir,outputdir,runid)	  
os.system(cmd)
  
  
#-------------------------------------------------------------



