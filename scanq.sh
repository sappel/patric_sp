#q=(4.01 4.03 4.05)
#for i in "${q[@]}"; do

sed -ie 's/Q_hor = .*/Q_hor = 4.48;/' mad/sis18_inj.mad;
sed -ie 's\subdir = ".*"  # subdirectory for output\subdir = "/q48_sc"  # subdirectory for output\' python/patric_2d_madx.py
python/patric_2d_madx.py
#done