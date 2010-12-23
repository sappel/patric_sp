# Call make with debug option and attach gdb to consecutively launched process
path=sim_inj/tmp  #  destination path

make MAKEOPTIONS=debug

#start simulation and gdb
python/patric_2d_madx.py  # should continue only if successful here, but how?
while [ ! -r ../$path/out.1.dat ] || [ $(wc -l < ../$path/out.1.dat) -eq 0 ]; do sleep 1; done
pid=$(perl -n -e "m%PID (\d+) % && print \"\$1\n\"" ../$path/out.1.dat | head -n 1)
gdb -x start.gdb patric $pid

