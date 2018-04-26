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

#create a topology with HMR
tmpfl=`mktemp /tmp/prep-parmed-XXXX`
cat <<EOF > $tmpfl
parm 148l_watio.prmtop
HMassRepartition
parmout 148l_watio_hmr.prmtop
go
EOF
parmed -i $tmpfl
rm $tmpfl

#Run equilibration
(cd Equil && ./equil.sh)

#Create an AMBER (DISANG) restraits file and update DU positions in the restart file.
python3 ../../FRETrest.py -t 148l_watio_hmr.prmtop -r Equil/26_md_nvt_red_pme_11.restrt -j t4l.fps.json --chi2 'C3 χ²' --fout prod_0001.f --restout prod_0000.restrt --force 50 --resoffset -1
ambpdb -p 148l_watio_hmr.prmtop -c prod_0000.restrt > prod_0000.pdb

#Run production
echo "Running production: "$( date --rfc-3339=s )
pmemd.cuda -O -i "prod_0001.in" -o "prod_0001.out" -p "148l_watio_hmr.prmtop" -c "prod_0000.restrt" -r "prod_0001.restrt" -x "prod_0001.mdcrd"
echo "Done: "$( date --rfc-3339=s )
ambpdb -p 148l_watio_hmr.prmtop -c prod_0001.restrt > prod_0001.pdb

# remove all the results
# rm leap.log 148l_PA.pdb prod_0001.pdb prod_0000.pdb prod_0000.restrt prod_0001.dump prod_0001.mdcrd prod_0001.out prod_0001.restrt 148l_watio.prmtop 148l_watio_hmr.prmtop 148l_upd.restrt prod_0001.f 148l_watio.inpcrd AV_*.pqr mdinfo
# rm Equil/logfile Equil/mdinfo Equil/*.out Equil/*.mdcrd Equil/*.restrt
