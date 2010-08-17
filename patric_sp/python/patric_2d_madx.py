#!/usr/bin/python
"""
A python script to launch PATRIC jobs 
"""

from math import *
from odict import * 
from socket import gethostname
import os

path='../../inj_sim/KV0/'  # path for output
runid=1            # run identification number (only out.dat)

#-----MPI options------------------

num_procs=2
machines=['lxir033']

#----------------------------------

patric_dict=OrderedDict()

patric_dict['NPIC']=10000  # 100000; particles per beamlett; SP
patric_dict['NX']=128
patric_dict['NY']=128
patric_dict['NZ']=256
patric_dict['NZ_bunch']=64
patric_dict['cells']=20
patric_dict['lossTol']=0.8  # tolerable relative losses; SP
patric_dict['e_kin']=11.4
patric_dict['Z']=28.0
patric_dict['A']=238.0
patric_dict['current']=0.015
patric_dict['piperadius']=0.1
patric_dict['coll_halfgap']=0.07  # radius of collimator representing the septum, normally equal to x_septum
patric_dict['circum']=216.0
patric_dict['gamma_t']=4.79
patric_dict['CF_advance_h']=97.3*pi/180.0
patric_dict['CF_advance_v']=97.3*pi/180.0 # 129.3*pi/180.0
patric_dict['CF_R']=0.0
patric_dict['CF_length']=18.0
patric_dict['NCF']=16
patric_dict['koct']=8.0
patric_dict['dQxm']=-0.0054*1.0
patric_dict['dQym']=-0.0054*1.0
patric_dict['dqx_detune']=-0.0025/6.0  # rms
patric_dict['dqy_detune']=-0.0025/6.0  # rms
patric_dict['pic_subset']=1  # 10000
patric_dict['init_pic_xy']=0  # 0 (WB), 1 (KV), 2 (SG), 3 (GS)
patric_dict['init_pic_z']=2  # 0 (elliptic, coast), 1 (elliptic bunch), 2 (Gauss, coast)  # was 4 !
patric_dict['momentum_spread']=1.0e-4  # rms
patric_dict['rms_emittance_x0']=1.3  # rms, mm mrad
patric_dict['rms_emittance_y0']=4.1  # rms, mm mrad
patric_dict['mismatch_x']=1.0  # initial beam size does not match beta function (1 = match)
patric_dict['mismatch_y']=1.0
patric_dict['x_septum']=0.07  # distance of septum from nominal orbit
patric_dict['offcenter']=0.0  # !!! displacement of barycenter from ideal injection [m]
patric_dict['inj_angle']=6.e-3  # injection angle in rad; added by SP
patric_dict['max_inj']=1  # maximal number of injections; added by SP
patric_dict['inj_phase']=0.15  # phase of injected beamletts in xx'-space; added by SP
patric_dict['sept_return']=5  # number of revolution till return of beamlett to septum, should be close to 1/Q_f [notes p. 139]; added by SP
patric_dict['bunchfactor']=1.0
patric_dict['dqci']=0.0
patric_dict['dqcr']=-0.0
patric_dict['Rs']=0.0 # 8.5e6
patric_dict['nres']=10.0
patric_dict['Qs']=2.0
patric_dict['leit']=1.0e5
patric_dict['Zimage']=-7.0e7 # -3.0e7/2.0
patric_dict['madx_input_file']=1  # 0 (no, i.e. use CF), 1 (yes)
patric_dict['space_charge']=0   #   0 (off), 1 (self-consistent), 2 (linear), 3 (nonlinear)
patric_dict['imp_kick']=0       #   0 (off), 1 (on)
patric_dict['sliced']=0
patric_dict['cavity']=0   # 0 (off), 1 (rf), 2 (barrier)
patric_dict['octupole_kick']=0
patric_dict['ampdetun_kick']=0  # commented out in Main.cpp, as incompatibel with modifications
patric_dict['chroma']=0
patric_dict['bc_end']=1  #  0 (open) , 1 (periodic)
patric_dict['print_cell']=1  # with MADX input one cell is the whole lattice
patric_dict['footprint']=0
patric_dict['btf']=0
patric_dict['btf_harmonic']=0
patric_dict['ausgabe']=path

#-------------------run-------------------------------------

if(os.path.exists(path)):
    rm = raw_input('Delete all files in directory '+path+'? (n/y) ')
    if(rm == 'y'):
        os.system('rm -r '+path)
    else:
        assert "Target directory could not be cleared."
os.mkdir(path)

# create new PATRIC configuration file
outfile=open('%s/patric.cfg' % (path),'w')
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
   cmd="nohup mpirun2 -np %d -machinefile machines.tmp ../bin/patric %s/patric > %s/out.%s.dat &" %(num_procs, path, path, runid)
else:
   cmd="nohup mpirun2 -nolocal -np %d -machinefile machines.tmp ../bin/patric %s/patric > %s/out.%s.dat &" %(num_procs, path, path, runid)
os.system(cmd)

#-------------------------------------------------------------