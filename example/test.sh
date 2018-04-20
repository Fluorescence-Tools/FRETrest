#!/bin/bash

python3 ../placeAVmp.py -p 148l_noH.pdb -o 148l_PA.pdb -j t4l.fps.json --chi2 'C3 χ²' --avprefix 'AV_'

#make AMBER input files
tmpfl=`mktemp /tmp/prep-leap-XXXX`
cat <<EOF > $tmpfl
source leaprc.protein.ff14SB
set default PBradii mbondi2
source leaprc.water.tip3p
loadamberparams DUM.frcmod
loadOff DUM.lib

PDB = loadpdb 148l_PA.pdb

addions PDB NA 0
addions PDB CL 0

solvateoct PDB TIP3PBOX 11.0
saveamberparm PDB 148l_watio.prmtop 148l_watio.inpcrd
quit
EOF

tleap -f $tmpfl
rm $tmpfl 

#rm leap.log 148l_PA.pdb 148l_watio.prmtop 148l_watio.inpcrd AV_*.pqr