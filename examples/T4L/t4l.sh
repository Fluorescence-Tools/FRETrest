#!/bin/bash

#append DUmmy atoms to a PDB. DU atoms represent AV mean positions.
python3 ../../placeAVmp.py -p 148l_noH.pdb -o 148l_PA.pdb -j t4l.fps.json --chi2 'C3 χ²' #--avprefix 'AV_'

#prepare AMBER input files
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

#TODO: add HMR

#Create an AMBER (DISANG) restraits file and update DU positions in the restart file.
python3 ../../FRETrest.py -t 148l_watio.prmtop -r Equil/11_md_nvt_ntr_pme.restrt -j t4l.fps.json --chi2 'C3 χ²' --fout 148l.f --restout fixed_11.restrt --force 50 --resoffset -1

#rm leap.log 148l_PA.pdb 148l_watio.prmtop 148l_upd.restrt 148l.f 148l_watio.inpcrd AV_*.pqr
