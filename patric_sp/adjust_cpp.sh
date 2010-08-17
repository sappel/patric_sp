# script to adjust spacing around '=', '==' and similar, before comments, after ',' and to delete spaces at line endings in cpp files
sed -e 's#\([^ =!<>+-]\)=#\1 =#g' -e 's/=\([^ =+-]\)/= \1/g' -e 's/ *\/\//  \/\//' -e 's/ *$//' -e 's/,\([^ \n]\)/, \1/g' -e 's/ ;/;/g' infile > outfile
