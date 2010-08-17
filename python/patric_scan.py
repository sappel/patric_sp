"""
A python script to launch PATRIC jobs 
"""
from math import *
from odict import * 
from socket import gethostname
import os

#outputdir='/d/bhs01/boine/track_sc_pb_n10Q2_fine/'
#outputdir='/d/bhs01/boine/track_sc_bb_n10Q2_fine_new/'
outputdir='/d/bhs01/boine/track_sc_dc_el_n10Q2_fine_new/'

#-----MPI options------------------

num_procs=2
#machines=['lxir002','lxir002','lxir003','lxir003','lxir004','lxir004','lxir005','lxir005']
#machines=['lxp025','lxp025','lxp026','lxp026','lxp027','lxp027','lxp028','lxp028']
machines=['lxir005','lxir005','lxir005','lxir005','lxir004','lxir007','lxir007','lxir007']

#----------------------------------

bunchf=1.0 # 0.8  
scalingf=10.0 

patric_dict=OrderedDict()

patric_dict['NPIC']=120000   # 120000
patric_dict['NX']=128
patric_dict['NY']=128
patric_dict['NZ']=128   # 128
patric_dict['NZ_bunch']=128
patric_dict['cells']=20480 # 16384 # 20480    
patric_dict['e_kin']=200.0
patric_dict['Z']=28.0
patric_dict['A']=238.0
patric_dict['current']=scalingf*0.7 # *0.666666    
patric_dict['piperadius']=0.04
patric_dict['circum']=216.0   # 54.0
patric_dict['gamma_t']=5.4
patric_dict['CF_advance_h']=97.3*pi/180.0
patric_dict['CF_advance_v']=97.3*pi/180.0
patric_dict['CF_R']=0.0
patric_dict['CF_length']=18.0
patric_dict['NCF']=16
patric_dict['koct']=8.0
patric_dict['dQxm']=-scalingf*0.0062/bunchf       #  -0.0135*2.0/0.8
patric_dict['dQym']=-scalingf*0.0062/bunchf       # -0.0135*2.0/0.8
patric_dict['dqx_detune']=0.04
patric_dict['dqy_detune']=0.04
patric_dict['pic_subset']=20000
patric_dict['init_pic_xy']=1
patric_dict['init_pic_z']=0  # 1 (bunch), 0 (coast), 2 (Gauss, coast)
patric_dict['momentum_spread']=scalingf*3.0e-4/bunchf  #  2.0*0.0005/0.8
patric_dict['rms_emittance_x0']=3.75
patric_dict['rms_emittance_y0']=3.75
patric_dict['mismatch_x']=1.0
patric_dict['mismatch_y']=1.0
patric_dict['offcenter']=0.001
patric_dict['bunchfactor']=bunchf
patric_dict['dqci']=0.0
patric_dict['dqcr']=0.0
patric_dict['Rs']=8.5e6
patric_dict['nres']=10.0
patric_dict['Qs']=2.0  
patric_dict['Zimage']=-1.0e7
patric_dict['madx_input_file']=0
patric_dict['space_charge']=2  
patric_dict['imp_kick']=1
patric_dict['sliced']=1
patric_dict['cavity']=0   # 0 (coast), 1 (rf), 2 (bb)
patric_dict['octupole_kick']=0
patric_dict['ampdetun_kick']=0
patric_dict['chroma']=0
patric_dict['bc_end']=1    # bunch: 0 , periodic: 1
patric_dict['print_cell']=120
patric_dict['footprint']=0
patric_dict['btf']=0
patric_dict['btf_harmonic']=0
patric_dict['ausgabe']='/d/bhs01/boine/track/run1'  

#-------------------loop-------------------------------------

# delete all old job output directories:	
#os.system("rm -rf %srun*" %(outputdir))

# machine counter	
machine_id=0

for jobid in range(1,9):   # 1-9
  # delete old directory (if exists):
  os.system("rm -rf %srun%d" %(outputdir,jobid))
  # create new job output directory:	
  os.system("mkdir %srun%d" %(outputdir,jobid))
  # create new PATRIC configuration file	
  outfile=open('%srun%d/patric.cfg' % (outputdir,jobid),'w')
  # change PATRIC input parameters: 
  patric_dict['ausgabe']='%srun%d' %(outputdir,jobid)  
  patric_dict['Rs']=machine_id*0.25e6 # 2.5e6      
  patric_dict['Zimage']=-3.0e7
  # write PATRIC configuration file in output dir:
  for x,y in patric_dict.iteritems():
	outfile.write('%s:  %s\n' %(x,y))
  outfile.close()
  # mpi machine configuation file:
  machinefile=open('tmp/machines.%d.tmp' % (jobid),'w')
  machinefile.write('%s:4\n' %(machines[machine_id]))
  machinefile.close()
  # start mpi job:
  thishost=gethostname()
  if thishost == machines[machine_id]:  
    cmd="nohup mpirun -np %d -machinefile tmp/machines.%d.tmp ../bin/track %srun%d/patric > tmp/out.%d.dat &" %(num_procs,jobid,outputdir,jobid,jobid)
  else:
    cmd="nohup mpirun -nolocal -np %d -machinefile tmp/machines.%d.tmp ../bin/track %srun%d/patric > tmp/out.%d.dat &" %(num_procs,jobid,outputdir,jobid,jobid)	  
  os.system(cmd)
  # next machine
  machine_id+=1 
  
#-------------------------------------------------------------



