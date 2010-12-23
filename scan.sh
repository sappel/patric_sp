q=(1 2 3 4 5 6 7 8)

j=0
while [ $j -lt 8 ]; do q[j]=$[${q[j]}+16]; j=$[j+1]; done

## using sed
#sed -ie "s/Q_hor = .*/Q_hor = 4.$i;/" mad/sis18_inj.mad;
#sed -i  \
#  -e 's\subdir = ".*"  \subdir = "/qKV19"  \'  \
#  python/patric_2d_madx.py
#python/patric_2d_madx.py


## using perl
for i in "${q[@]}"; do
    perl -pi -e "s/Q_hor = .*/Q_hor = 4.$(printf %02d $i);/" mad/sis18_inj.mad
    perl -pi -e "s\subdir = .*  \subdir = '/qKV$(printf %02d $i)_sc'  \\"   python/patric_2d_madx.py
    python/patric_2d_madx.py;
    sleep 0.1s
done


## fragment for later use
#  -e "s|offcenter']=.*|offcenter']=-0.000|" \
