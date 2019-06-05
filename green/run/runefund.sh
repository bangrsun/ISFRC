#!/bin/sh
#----------------------------------------------------------------------------
# Written by lzp 2015/11/12 for New EFUND code running as old style
#----------------------------------------------------------------------------

MHDIN=EXL_enn
TBLDIR=../../green_table
EXEDIR=.

# Removing existed data files
if [ -e mhdout.dat ]; then
	rm mhdout.dat
fi

# run efundud
echo "cp ${MHDIN} mhdin.dat"
cp ${MHDIN} mhdin.dat
echo "${EXEDIR}/efundud"
${EXEDIR}/efundud

# generate dprobe.dat
ex mhdin.dat <<EOF
2,12 d
wq! dprobe.dat
EOF
# move output of efundud to table directory
echo "mv *.ddd dprobe.dat ${TBLDIR}"
mv *.ddd dprobe.dat ${TBLDIR}

echo "The Green Tables have moved to green_table directory!"
