# call MADX if necessary
if [  mad/sis18_inj.mad -nt mad/twiss_inj.txt  -o  mad/sis18_inj.mad -nt mad/sectormap_inj.txt ];
then 
    cd mad/
    madx<sis18_inj.mad
    cd ..;
fi

# clean output directory
files=(dipole_kick_x.dat dipole_x.dat dipole_y.dat Ex.dat Ey.dat idl.dat out.1.dat patric.cfg patric.dat pics_0.dat pics_1.dat rho_xy.dat rho_z.dat xsys.dat xxs.dat yys.dat zx.dat)
for i in ${files[*]}; do
    if [ -e ../output/"$i" ] 
	then rm ../output/"$i"
    fi;
done

# call make with debug option
make MAKEOPTIONS=debug

#start simulation and gdb
cd ./python/
python patric_2d_madx.py
cd ../
while [ ! -r ../output/out.1.dat ] || [ $(wc -l < ../output/out.1.dat) -eq 0 ]; do sleep 1; done
pid=$(perl -n -e "m%PID (\d+) % && print \"\$1\n\"" ../output/out.1.dat | head -n 1)
gdb -x start.gdb patric $pid
